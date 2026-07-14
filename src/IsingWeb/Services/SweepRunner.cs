using System.Diagnostics;
using System.Text;

using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace IsingWeb.Services;

public sealed record SweepPoint(
    double Temperature,
    double MagnetisationAbsolute,
    double MagnetisationAbsoluteSigma,
    double Energy,
    double EnergySigma,
    double Susceptibility,
    double SusceptibilitySigma,
    double SpecificHeat,
    double SpecificHeatSigma,
    double CorrelationLengthOverL,
    double CorrelationLengthOverLSigma);

public enum SweepQuality
{
    Quick,
    Standard,
    Thorough
}

/// <summary>
///     Runs an observables-vs-temperature sweep in browser-friendly chunks,
///     yielding to the event loop so the UI stays responsive.
/// </summary>
public sealed class SweepRunner
{
    private CancellationTokenSource? _cancellation;

    public bool IsRunning { get; private set; }

    public double Progress { get; private set; }

    public string Status { get; private set; } = "";

    public List<SweepPoint> Points { get; } = new();

    public event Action? Updated;

    public void Cancel() => _cancellation?.Cancel();

    public async Task RunAsync(
        int latticeLength,
        double minTemperature,
        double maxTemperature,
        int temperatureCount,
        SpinUpdateMethod method,
        SweepQuality quality)
    {
        if (IsRunning)
        {
            return;
        }

        _cancellation = new CancellationTokenSource();
        var token = _cancellation.Token;
        IsRunning = true;
        Points.Clear();
        Progress = 0.0;

        var (thermalisationSweeps, sweepsPerMeasurement, measurementsPerBlock, blocks) = quality switch
        {
            SweepQuality.Quick => (150, 2, 25, 2),
            SweepQuality.Thorough => (800, 8, 60, 4),
            _ => (400, 4, 40, 3)
        };

        var totalSpins = latticeLength * latticeLength;
        var temperatures = Enumerable.Range(0, temperatureCount)
                                     .Select(i => minTemperature
                                                + (maxTemperature - minTemperature) * i
                                                / Math.Max(1, temperatureCount - 1))
                                     .ToList();

        try
        {
            for (var index = 0; index < temperatures.Count; index++)
            {
                var temperature = temperatures[index];
                var beta = 1.0 / temperature;
                Status = $"T = {temperature:0.000} ({index + 1}/{temperatures.Count}): thermalising…";
                NotifyProgress(index, temperatures.Count, phase: 0.0);

                var initialSpins = new List<int>(totalSpins);
                for (var i = 0; i < totalSpins; i++)
                {
                    initialSpins.Add(Random.Shared.Next(2) is 0 ? 1 : -1);
                }

                var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength, initialSpins);

                await RunChunkedAsync(
                    simulator,
                    beta,
                    (long)thermalisationSweeps * totalSpins,
                    method,
                    token);

                simulator.FlushWolffQueue(beta, SimulationRunner.FerromagneticJ, h: 0.0);

                Status = $"T = {temperature:0.000} ({index + 1}/{temperatures.Count}): measuring…";
                NotifyProgress(index, temperatures.Count, phase: 0.3);

                var blockMagnetisationAbs = new List<double>();
                var blockEnergy = new List<double>();
                var blockSusceptibility = new List<double>();
                var blockSpecificHeat = new List<double>();
                var blockXiOverL = new List<double>();

                for (var block = 0; block < blocks; block++)
                {
                    simulator.ResetObservables();
                    BasicObservables? observables = null;

                    // Chunk at measurement granularity: without a reset the
                    // simulator's accumulators keep running, so the last call
                    // returns the running average over the whole block.
                    var chunkStopwatch = Stopwatch.StartNew();
                    for (var m = 0; m < measurementsPerBlock; m++)
                    {
                        token.ThrowIfCancellationRequested();
                        observables = simulator.RunMonteCarloWithObservablesComputation(
                            beta,
                            SimulationRunner.FerromagneticJ,
                            h: 0.0,
                            iterationsNeededForSingleChiXiMeasurement: sweepsPerMeasurement * totalSpins,
                            measurementsCountForChiXiExpectationValue: 1,
                            method);

                        if (chunkStopwatch.Elapsed.TotalMilliseconds > 30)
                        {
                            NotifyProgress(
                                index,
                                temperatures.Count,
                                phase: 0.3
                                     + 0.7 * (block + (double)(m + 1) / measurementsPerBlock) / blocks);
                            await Task.Delay(1, token);
                            chunkStopwatch.Restart();
                        }
                    }

                    if (observables is null)
                    {
                        continue;
                    }

                    blockMagnetisationAbs.Add(observables.MagnetisationAbsolute);
                    blockEnergy.Add(observables.Energy);
                    blockSusceptibility.Add(observables.Susceptibility);
                    blockSpecificHeat.Add(observables.SpecificHeat);
                    if (!double.IsNaN(observables.RenormalisedCorrelationLength))
                    {
                        blockXiOverL.Add(observables.RenormalisedCorrelationLength);
                    }
                }

                Points.Add(new SweepPoint(
                    temperature,
                    Mean(blockMagnetisationAbs),
                    Sigma(blockMagnetisationAbs),
                    Mean(blockEnergy),
                    Sigma(blockEnergy),
                    Mean(blockSusceptibility),
                    Sigma(blockSusceptibility),
                    Mean(blockSpecificHeat),
                    Sigma(blockSpecificHeat),
                    Mean(blockXiOverL),
                    Sigma(blockXiOverL)));

                NotifyProgress(index + 1, temperatures.Count, phase: 0.0);
            }

            Status = "Done.";
            Progress = 1.0;
        }
        catch (OperationCanceledException)
        {
            Status = "Cancelled.";
        }
        finally
        {
            IsRunning = false;
            _cancellation.Dispose();
            _cancellation = null;
            Updated?.Invoke();
        }
    }

    public string ToCsv()
    {
        var builder = new StringBuilder();
        builder.AppendLine(
            "T,abs_m,abs_m_sigma,energy_per_site,energy_sigma,susceptibility,susceptibility_sigma,"
          + "specific_heat,specific_heat_sigma,xi_over_L,xi_over_L_sigma");

        foreach (var p in Points)
        {
            builder.AppendLine(
                $"{p.Temperature},{p.MagnetisationAbsolute},{p.MagnetisationAbsoluteSigma},"
              + $"{p.Energy},{p.EnergySigma},{p.Susceptibility},{p.SusceptibilitySigma},"
              + $"{p.SpecificHeat},{p.SpecificHeatSigma},"
              + $"{p.CorrelationLengthOverL},{p.CorrelationLengthOverLSigma}");
        }

        return builder.ToString();
    }

    private async Task RunChunkedAsync(
        IsingMonteCarloSimulator simulator,
        double beta,
        long totalSteps,
        SpinUpdateMethod method,
        CancellationToken token)
    {
        var chunk = Math.Max(1024, simulator.TotalSpinsCount);
        long done = 0;
        var stopwatch = Stopwatch.StartNew();

        while (done < totalSteps)
        {
            token.ThrowIfCancellationRequested();
            var stepsNow = Math.Min(chunk, totalSteps - done);
            simulator.RunSteps(stepsNow, beta, SimulationRunner.FerromagneticJ, h: 0.0, method);
            done += stepsNow;

            if (stopwatch.Elapsed.TotalMilliseconds > 30)
            {
                await Task.Delay(1, token);
                stopwatch.Restart();
            }
        }
    }

    private void NotifyProgress(double completedTemperatures, int totalTemperatures, double phase)
    {
        Progress = Math.Clamp((completedTemperatures + phase) / totalTemperatures, 0.0, 1.0);
        Updated?.Invoke();
    }

    private static double Mean(IReadOnlyCollection<double> values) =>
        values.Count is 0 ? double.NaN : values.Average();

    private static double Sigma(IReadOnlyCollection<double> values)
    {
        if (values.Count < 2)
        {
            return double.NaN;
        }

        var mean = values.Average();
        return Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1.0));
    }
}
