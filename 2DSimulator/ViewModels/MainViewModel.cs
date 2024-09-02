using System;
using System.Collections.Generic;
using System.Linq;
using System.Timers;

using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;

using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace _2DSimulator.ViewModels;

public partial class MainViewModel : ViewModelBase
{
    private const int RandomSeed = 41;
    private const int FixedLatticeLength = 128;
    private const int StepsCount = 200;

    private readonly Timer timer = new(interval: 50) { Enabled = false };

    [ObservableProperty]
    [NotifyPropertyChangedFor(nameof(Title))]
    public double magnetisation;

    public string Title => $"2D Ising simulation {LatticeLength} x {LatticeLength}  //  M = {Magnetisation}";

    public string Greeting => "Hi!";

    public int LatticeLength { get; }

    public Spin[][] Spins { get; set; }

    [ObservableProperty]
    private string toggleText;

    [ObservableProperty]
    private double energy;

    private readonly int _totalSpinsCount;
    private readonly Random _random;
    private readonly Queue<int[]> _clusterQueue;

    private int _referenceSpin;

    private readonly double _temperature;
    private readonly double _j;
    private readonly double _h;
    private readonly double? _jY;
    private readonly SpinUpdateMethod _spinUpdateMethod;

    public MainViewModel()
        : this(FixedLatticeLength) { }

    public MainViewModel(int latticeLength)
    {
        LatticeLength = latticeLength;
        _totalSpinsCount = LatticeLength * LatticeLength;

        InitialSpinConfiguration = Fast2DIsingMonteCarloSimulator.Get2DArrayFrom1DArray(
            SpinConfigurationBuilder.InitialiseRandomLattice(
                _totalSpinsCount,
                initialSpinDownRatio: 0.25,
                RandomSeed));

        ToSpinsWithSpinObject();

        _j = -1.0;
        _jY = null;
        _h = -0.0;

        var totalMagnetisation = Convert.ToDouble(InitialSpinConfiguration.SelectMany(elements => elements).Sum());
        Magnetisation = totalMagnetisation
                      / _totalSpinsCount;

        Energy = 0;

        for (var y = 0; y < LatticeLength; y++)
        {
            for (var x = 0; x < LatticeLength; x++)
            {
                Energy += 0.5 * GetEnergyOfSite(x, y, _j, _h, _jY);
            }
        }

        _random = new Random(RandomSeed);
        _clusterQueue = new Queue<int[]>();

        _spinUpdateMethod = SpinUpdateMethod.Glauber;
        _temperature = 2.0;

        timer.Elapsed += (_, _) => Update(1.0 / _temperature, _spinUpdateMethod, StepsCount);
    }

    public int[][] InitialSpinConfiguration { get; set; }

    private void ToSpinsWithSpinObject()
    {
        Spins = new Spin[LatticeLength][];
        for (var row = 0; row < LatticeLength; row++)
        {
            Spins[row] = new Spin[LatticeLength];
            for (var column = 0; column < LatticeLength; column++)
            {
                Spins[row][column] = new Spin(row, column, InitialSpinConfiguration[row][column]);
            }
        }
    }

    [RelayCommand]
    private void Randomize()
    {
        var timerEnabled = timer.Enabled;
        timer.Enabled = false;

        InitialSpinConfiguration = Fast2DIsingMonteCarloSimulator.Get2DArrayFrom1DArray(
            SpinConfigurationBuilder.InitialiseRandomLattice(
                _totalSpinsCount,
                _random.NextDouble()));

        ToSpinsWithSpinObject();

        timer.Enabled = timerEnabled;
    }

    public void Update(double givenBeta, SpinUpdateMethod givenSpinUpdateMethod, int stepsCount)
    {
        for (var i = 0; i < stepsCount; i++)
        {
            if (givenSpinUpdateMethod is SpinUpdateMethod.Glauber)
            {
                var chosenX = _random.Next(minValue: 0, LatticeLength);
                var chosenY = _random.Next(minValue: 0, LatticeLength);

                if (FlipWithGlauber(_random.NextDouble(), chosenX, chosenY, givenBeta))
                {
                    FlipSpin(chosenX, chosenY, _j, _h, _jY);
                }
            }
            else
            {
                FlipSpinWithWolff(givenBeta);
            }
        }
    }

    public IEnumerable<Spin> All
    {
        get
        {
            var rows = LatticeLength;
            var columns = LatticeLength;

            for (var r = 0; r < rows; r++)
            {
                for (var c = 0; c < columns; c++)
                {
                    yield return Spins[r][c];
                }
            }
        }
    }

    [RelayCommand]
    public void ToggleRun(bool? state = null)
    {
        timer.Enabled = state ?? !timer.Enabled;
        ToggleText = timer.Enabled ? "⏸ Pause" : "▶ Run";
    }

    public void Dispose() => timer.Dispose();

    private bool FlipWithGlauber(
        double randomProbability,
        int x,
        int y,
        double beta) =>
        randomProbability <= 1.0 / (1.0 + Math.Exp(beta * GetDeltaEnergyOfSite(x, y, _j, _h, _jY)));

    private void FlipSpinWithWolff(double beta)
    {
        // If a cluster has not been chosen yet
        if (_clusterQueue.Count is 0)
        {
            var chosenX = _random.Next(minValue: 0, LatticeLength);
            var chosenY = _random.Next(minValue: 0, LatticeLength);
            _referenceSpin = Spins[chosenX][chosenY].SpinValue;
            _clusterQueue.Enqueue(new[] { chosenX, chosenY });

            // Flip reference spin of the cluster
            FlipSpin(chosenX, chosenY, _j, _h, _jY);
        }

        var siteIndexToConsider = _clusterQueue.Dequeue();
        var neighboursIndices = GetNearestNeighboursIndices(siteIndexToConsider[0], siteIndexToConsider[1]);

        foreach (var neighbour in neighboursIndices)
        {
            if (Spins[neighbour[0]][neighbour[1]].SpinValue != _referenceSpin)
            {
                continue;
            }

            var addToCluster = _random.NextDouble() <= 1.0 - Math.Exp(-2.0 * beta);
            if (addToCluster)
            {
                _clusterQueue.Enqueue(neighbour);
                FlipSpin(neighbour[0], neighbour[1], _j, _h, _jY);
            }
        }
    }

    public void FlipSpin(int x, int y, double j, double h, double? jY)
    {
        Energy += GetDeltaEnergyOfSite(x, y, j, h, jY);
        Magnetisation += Convert.ToDouble(-2 * Spins[x][y].SpinValue) / _totalSpinsCount;

        Spins[x][y].SpinValue *= -1;
    }

    public double GetDeltaEnergyOfSite(
        int x,
        int y,
        double j,
        double h,
        double? jY = null) =>
        -2.0 * GetEnergyOfSite(x, y, j, h, jY);

    public double GetEnergyOfSite(
        int x,
        int y,
        double j,
        double h,
        double? jY = null)
    {
        var spinValue = Spins[x][y].SpinValue;
        var neighbours = GetNearestNeighboursIndices(x, y);

        // Possible in 2D only
        if (jY is not null)
        {
            var xBonds = j
                       * neighbours.Take(count: 2)
                                   .Select(
                                       neighbourIndex =>
                                           spinValue
                                         * Spins[neighbourIndex[0]][neighbourIndex[1]].SpinValue)
                                   .Sum();
            var yBonds = jY
                       * neighbours.TakeLast(count: 2)
                                   .Select(
                                       neighbourIndex =>
                                           spinValue
                                         * Spins[neighbourIndex[0]][neighbourIndex[1]].SpinValue)
                                   .Sum();

            return xBonds + (double)yBonds - h * spinValue;
        }

        return j
             * neighbours
               .Select(
                   neighbourIndex => spinValue
                                   * Spins[neighbourIndex[0]][neighbourIndex[1]].SpinValue)
               .Sum()
             - h * spinValue;
    }

    private int[][] GetNearestNeighboursIndices(int x, int y) =>
    [
        [(x + LatticeLength - 1) % LatticeLength, y], // Left
        [(x + 1) % LatticeLength, y], // Right
        [x, (y + LatticeLength - 1) % LatticeLength], // Up
        [x, (y + 1) % LatticeLength] // Down
    ];
}
