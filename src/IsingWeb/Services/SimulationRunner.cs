using System.Diagnostics;

using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace IsingWeb.Services;

/// <summary>
///     Owns the interactive simulation: parameters, the frame-budgeted stepping
///     loop, and the live M(t)/E(t) series. UI components read and mutate it;
///     changes take effect on the next animation-frame tick.
/// </summary>
public sealed class SimulationRunner
{
    public const double FerromagneticJ = -1.0;

    private const int SeriesCapacity = 400;

    private readonly Stopwatch _rateStopwatch = Stopwatch.StartNew();

    private long _stepsSinceRateSample;
    private double _lastSampleSweep;
    private double _field;

    // Fractional sweeps carried over between 1D frames so that speeds below
    // one sweep per frame skip frames instead of flooring to one row.
    private double _sweepCarry;
    private SpinUpdateMethod _method = SpinUpdateMethod.Metropolis;

    public SimulationRunner() => Configure(dimension: 2, latticeLength: 256, hotStart: true);

    public IsingMonteCarloSimulator Simulator { get; private set; } = null!;

    public int Dimension { get; private set; }

    public int LatticeLength { get; private set; }

    public bool HotStart { get; private set; }

    public double Temperature { get; set; } = ExactIsing2D.CriticalTemperature;

    public bool IsRunning { get; set; }

    /// <summary>Target sweeps per animation frame; the ~10 ms frame budget caps it.</summary>
    public double SweepsPerFrame { get; set; } = 2.0;

    public double CompletedSweeps { get; private set; }

    public double SweepsPerSecond { get; private set; }

    public LatticeColourMode ColourMode { get; set; } = LatticeColourMode.Spins;

    public bool HighlightWolffCluster { get; set; } = true;

    /// <summary>Fires after each completed sweep while in 1D mode (spacetime row).</summary>
    public Action? SweepCompletedIn1D { get; set; }

    public IReadOnlyList<double> SweepAxis => _sweepAxis;

    public IReadOnlyList<double> MagnetisationSeries => _magnetisationSeries;

    public IReadOnlyList<double> EnergySeries => _energySeries;

    private readonly List<double> _sweepAxis = new();
    private readonly List<double> _magnetisationSeries = new();
    private readonly List<double> _energySeries = new();

    public double Field
    {
        get => _field;
        set
        {
            // The library rejects Wolff with a field; the UI also disables it.
            _field = Method is SpinUpdateMethod.Wolff ? 0.0 : Math.Clamp(value, -2.0, 2.0);
            Simulator.Hamiltonian.InvalidateTotalEnergy();
        }
    }

    public SpinUpdateMethod Method
    {
        get => _method;
        set
        {
            if (value == _method)
            {
                return;
            }

            // Never abandon a half-grown cluster when leaving Wolff.
            Simulator.FlushWolffQueue(1.0 / Temperature, FerromagneticJ, _field);
            _method = value;

            if (_method is SpinUpdateMethod.Wolff)
            {
                _field = 0.0;
                Simulator.Hamiltonian.InvalidateTotalEnergy();
            }
        }
    }

    public double Magnetisation => Simulator.Hamiltonian.GetAverageMagnetisation(FerromagneticJ, _field);

    public double EnergyPerSite => Simulator.Hamiltonian.GetAverageEnergy(FerromagneticJ, _field);

    public void Configure(int dimension, int latticeLength, bool hotStart)
    {
        Dimension = dimension;
        LatticeLength = latticeLength;
        HotStart = hotStart;

        var totalSpins = dimension is 1 ? latticeLength : latticeLength * latticeLength;
        List<int>? initialSpins = null;
        if (hotStart)
        {
            initialSpins = new List<int>(totalSpins);
            for (var i = 0; i < totalSpins; i++)
            {
                initialSpins.Add(Random.Shared.Next(2) is 0 ? 1 : -1);
            }
        }

        Simulator = new IsingMonteCarloSimulator(dimension, latticeLength, initialSpins);
        CompletedSweeps = 0.0;
        _lastSampleSweep = 0.0;
        _stepsSinceRateSample = 0;
        _sweepCarry = 0.0;
        _sweepAxis.Clear();
        _magnetisationSeries.Clear();
        _energySeries.Clear();
    }

    public void Reset() => Configure(Dimension, LatticeLength, HotStart);

    /// <summary>
    ///     Advances the simulation for at most <paramref name="budgetMilliseconds" />
    ///     of wall time, aiming for <see cref="SweepsPerFrame" /> sweeps.
    /// </summary>
    public void Tick(double budgetMilliseconds)
    {
        if (!IsRunning)
        {
            return;
        }

        var beta = 1.0 / Temperature;
        var totalSpins = Simulator.TotalSpinsCount;
        var stopwatch = Stopwatch.StartNew();

        if (Dimension is 1)
        {
            // One spacetime row per completed sweep. Sub-1 rates accumulate in
            // _sweepCarry so slow speeds draw a row only every few frames.
            _sweepCarry += SweepsPerFrame;
            while (_sweepCarry >= 1.0)
            {
                Simulator.RunSteps(totalSpins, beta, FerromagneticJ, _field, _method);
                CompletedSweeps += 1.0;
                _stepsSinceRateSample += totalSpins;
                SweepCompletedIn1D?.Invoke();
                _sweepCarry -= 1.0;

                if (stopwatch.Elapsed.TotalMilliseconds > budgetMilliseconds)
                {
                    break;
                }
            }
        }
        else
        {
            var targetSteps = (long)(SweepsPerFrame * totalSpins);
            var batch = Math.Clamp(totalSpins / 8, 256, 65_536);
            long done = 0;
            while (done < targetSteps)
            {
                var stepsNow = Math.Min(batch, targetSteps - done);
                Simulator.RunSteps(stepsNow, beta, FerromagneticJ, _field, _method);
                done += stepsNow;

                if (stopwatch.Elapsed.TotalMilliseconds > budgetMilliseconds)
                {
                    break;
                }
            }

            CompletedSweeps += (double)done / totalSpins;
            _stepsSinceRateSample += done;
        }

        UpdateRate();
        SampleSeries();
    }

    public void StepOneSweep()
    {
        var beta = 1.0 / Temperature;
        Simulator.RunSteps(Simulator.TotalSpinsCount, beta, FerromagneticJ, _field, _method);
        Simulator.FlushWolffQueue(beta, FerromagneticJ, _field);
        CompletedSweeps += 1.0;
        if (Dimension is 1)
        {
            SweepCompletedIn1D?.Invoke();
        }

        SampleSeries(force: true);
    }

    public void Pause()
    {
        IsRunning = false;
        Simulator.FlushWolffQueue(1.0 / Temperature, FerromagneticJ, _field);
    }

    private void UpdateRate()
    {
        var elapsed = _rateStopwatch.Elapsed.TotalSeconds;
        if (elapsed < 0.5)
        {
            return;
        }

        SweepsPerSecond = _stepsSinceRateSample / (double)Simulator.TotalSpinsCount / elapsed;
        _stepsSinceRateSample = 0;
        _rateStopwatch.Restart();
    }

    private void SampleSeries(bool force = false)
    {
        if (!force && CompletedSweeps - _lastSampleSweep < 1.0)
        {
            return;
        }

        _lastSampleSweep = CompletedSweeps;
        _sweepAxis.Add(CompletedSweeps);
        _magnetisationSeries.Add(Magnetisation);
        _energySeries.Add(EnergyPerSite);

        if (_sweepAxis.Count > SeriesCapacity)
        {
            _sweepAxis.RemoveAt(0);
            _magnetisationSeries.RemoveAt(0);
            _energySeries.RemoveAt(0);
        }
    }
}
