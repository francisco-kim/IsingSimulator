using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;

namespace IsingMonteCarlo.Representations;

public class Fast2DIsingMonteCarloSimulator
{
    private const int Hundred = 100;

    private readonly int _latticeLength;
    private readonly int _totalSpinsCount;
    private readonly Random _random;
    private readonly Queue<int[]> _clusterQueue;

    private int _referenceSpin;

    public Fast2DIsingMonteCarloSimulator(IEnumerable<int> spinConfiguration, int? randomSeed = null)
    {
        _totalSpinsCount = spinConfiguration.Count();
        _latticeLength = Convert.ToInt32(Math.Sqrt(spinConfiguration.Count()));

        if (_totalSpinsCount != _latticeLength * _latticeLength)
        {
            throw new ArgumentNullException(nameof(spinConfiguration));
        }

        SpinConfiguration = Get2DArrayFrom1DArray(spinConfiguration);

        _random = randomSeed is not null ? new Random((int)randomSeed) : new Random();
        _clusterQueue = new Queue<int[]>();
    }

    public Fast2DIsingMonteCarloSimulator(int[][] spinConfiguration, int? randomSeed = null)
    {
        _totalSpinsCount = spinConfiguration.Count();
        _latticeLength = Convert.ToInt32(Math.Sqrt(spinConfiguration.Count()));

        if (_totalSpinsCount != _latticeLength * _latticeLength)
        {
            throw new ArgumentNullException(nameof(spinConfiguration));
        }

        SpinConfiguration = spinConfiguration ?? throw new ArgumentNullException(nameof(spinConfiguration));

        _random = randomSeed is not null ? new Random((int)randomSeed) : new Random();
        _clusterQueue = new Queue<int[]>();
    }

    public int[][] SpinConfiguration { get; }

    public void RunMonteCarlo(
        double beta,
        long? iterationLimit,
        SpinUpdateMethod spinUpdateMethod,
        bool saveLattice)
    {
        if (iterationLimit is null)
        {
            while (true)
            {
                flipSpinWithChosenDynamics(beta, spinUpdateMethod);
            }
        }

        long iterationCount = 0;
        while (iterationCount < iterationLimit)
        {
            flipSpinWithChosenDynamics(beta, spinUpdateMethod);

            ++iterationCount;
        }

        if (spinUpdateMethod is SpinUpdateMethod.Wolff && _clusterQueue.Count != 0)
        {
            while (_clusterQueue.Count > 0)
            {
                FlipSpinWithWolff(beta);
            }
        }

        if (saveLattice)
        {
            SaveSpinConfiguration(1.0 / beta, (long)iterationLimit);
        }

        void flipSpinWithChosenDynamics(double givenBeta, SpinUpdateMethod givenSpinUpdateMethod)
        {
            if (givenSpinUpdateMethod is SpinUpdateMethod.Glauber)
            {
                var chosenX = _random.Next(minValue: 0, _latticeLength);
                var chosenY = _random.Next(minValue: 0, _latticeLength);

                if (FlipWithGlauber(randomProbability: _random.NextDouble(), chosenX, chosenY, givenBeta))
                {
                    FlipSpin(chosenX, chosenY);
                }
            }
            else
            {
                FlipSpinWithWolff(givenBeta);
            }
        }
    }

    public void RunVerboseMonteCarlo(
        double beta,
        long? iterationLimit,
        SpinUpdateMethod spinUpdateMethod,
        bool saveLattice)
    {
        if (iterationLimit is null)
        {
            while (true)
            {
                flipSpinWithChosenDynamics(beta, spinUpdateMethod);
            }
        }

        long iterationCount = 0;
        while (iterationCount < iterationLimit)
        {
            flipSpinWithChosenDynamics(beta, spinUpdateMethod);

            ++iterationCount;

            if (iterationCount % _totalSpinsCount == 0)
            {
                Console.Write($"{iterationCount / _totalSpinsCount} MC sweep completed.");
                Console.Write('\r');
            }
        }

        if (spinUpdateMethod is SpinUpdateMethod.Wolff && _clusterQueue.Count != 0)
        {
            while (_clusterQueue.Count > 0)
            {
                FlipSpinWithWolff(beta);
            }
        }

        if (saveLattice)
        {
            SaveSpinConfiguration(1.0 / beta, (long)iterationLimit);
        }

        void flipSpinWithChosenDynamics(double beta, SpinUpdateMethod spinUpdateMethod)
        {
            if (spinUpdateMethod is SpinUpdateMethod.Glauber)
            {
                var chosenX = _random.Next(minValue: 0, _latticeLength);
                var chosenY = _random.Next(minValue: 0, _latticeLength);

                if (FlipWithGlauber(randomProbability: _random.NextDouble(), chosenX, chosenY, beta))
                {
                    FlipSpin(chosenX, chosenY);
                }
            }
            else
            {
                FlipSpinWithWolff(beta);
            }
        }
    }

    public static List<int> ThermaliseLargeLattice(
        string? filename,
        int latticeLength,
        double temperature,
        long thermalisationStepsIn100MCSweepUnit = 1,
        SpinUpdateMethod spinUpdateMethod = SpinUpdateMethod.Wolff,
        bool resetIterationCountDuringSave = false,
        int? randomSeed = null,
        bool verbose = true)
    {
        var (initialSpinConfiguration, _, previousIterationCount) =
            SpinConfigurationBuilder.InitialiseLattice(
                filename,
                dimension: 2,
                latticeLength,
                initialSpinDownRatio: 0.25,
                randomSeed: randomSeed);

        var simulation = new Fast2DIsingMonteCarloSimulator(initialSpinConfiguration, randomSeed);

        if (verbose)
        {
            simulation.RunVerboseMonteCarlo(
                1.0 / temperature,
                thermalisationStepsIn100MCSweepUnit * Hundred * latticeLength * latticeLength,
                spinUpdateMethod,
                saveLattice: false);
        }
        else
        {
            simulation.RunMonteCarlo(
                1.0 / temperature,
                thermalisationStepsIn100MCSweepUnit * Hundred * latticeLength * latticeLength,
                spinUpdateMethod,
                saveLattice: false);
        }

        var iterationSteps = resetIterationCountDuringSave
                                 ? thermalisationStepsIn100MCSweepUnit
                                 : previousIterationCount
                                 + thermalisationStepsIn100MCSweepUnit;
        FileHelpers.SaveSpinConfiguration(
            simulation.GetSpinConfiguration(),
            temperature,
            iterationSteps,
            isInByte: true);

        Console.WriteLine(
            $"Thermalisation done.\n");

        return simulation.GetSpinConfiguration();
    }

    private bool FlipWithGlauber(
        double randomProbability,
        int x,
        int y,
        double beta) =>
        randomProbability <= 1.0 / (1.0 + Math.Exp(beta * GetDeltaEnergyOfSite(x, y)));

    private void FlipSpinWithWolff(double beta)
    {
        // If a cluster has not been chosen yet
        if (_clusterQueue.Count is 0)
        {
            var chosenX = _random.Next(minValue: 0, _latticeLength);
            var chosenY = _random.Next(minValue: 0, _latticeLength);
            _referenceSpin = SpinConfiguration[chosenX][chosenY];
            _clusterQueue.Enqueue(new[] { chosenX, chosenY });

            // Flip reference spin of the cluster
            FlipSpin(chosenX, chosenY);
        }

        var siteIndexToConsider = _clusterQueue.Dequeue();
        var neighboursIndices = GetNearestNeighboursIndices(siteIndexToConsider[0], siteIndexToConsider[1]);

        foreach (var neighbour in neighboursIndices)
        {
            if (SpinConfiguration[neighbour[0]][neighbour[1]] != _referenceSpin)
            {
                continue;
            }

            var addToCluster = _random.NextDouble() <= 1.0 - Math.Exp(-2.0 * beta);
            if (addToCluster)
            {
                _clusterQueue.Enqueue(neighbour);
                FlipSpin(neighbour[0], neighbour[1]);
            }
        }
    }

    public void FlipSpin(int x, int y) => SpinConfiguration[x][y] *= -1;

    public double GetDeltaEnergyOfSite(int x, int y) => -2.0 * GetEnergyOfSite(x, y);

    public double GetEnergyOfSite(int x, int y)
    {
        var spinValue = SpinConfiguration[x][y];

        return -1 * GetNearestNeighboursIndices(x, y).Select(neighbourIndex => SpinConfiguration[neighbourIndex[0]][neighbourIndex[1]]
                                                                        * spinValue)
                                              .Sum();
    }

    public List<int> GetSpinConfiguration()
    {
        var lattice = new List<int>(_totalSpinsCount);
        for (var y = 0; y < _latticeLength; ++y)
        {
            for (var x = 0; x < _latticeLength; ++x)
            {
                lattice.Add(SpinConfiguration[x][y]);
            }
        }

        return lattice;
    }

    private int[][] GetNearestNeighboursIndices(int x, int y) => new int[][]
        {
            new[] { (x + _latticeLength - 1) % _latticeLength, y }, // Left
            new[] { (x + 1) % _latticeLength, y }, // Right
            new[] { x, (y + _latticeLength - 1) % _latticeLength }, // Up
            new[] { x, (y + 1) % _latticeLength } // Down
        };

    private int[][] Get2DArrayFrom1DArray(IEnumerable<int> spinConfiguration)
    {
        var spins = spinConfiguration?.ToList() ?? throw new ArgumentNullException(nameof(spinConfiguration));

        var lattice = new int[_latticeLength][];
        for (var i = 0; i < _latticeLength; ++i)
        {
            lattice[i] = new int[_latticeLength];
        }

        for (var y = 0; y < _latticeLength; ++y)
        {
            for (var x = 0; x < _latticeLength; ++x)
            {
                lattice[x][y] = spins[x + y * _latticeLength];
            }
        }

        return lattice;
    }

    private void SaveSpinConfiguration(double temperature, long iterationCount)
    {
        var latticeConfiguration = GetSpinConfiguration();

        FileHelpers.SaveSpinConfiguration(
            latticeConfiguration,
            temperature,
            iterationCount,
            isInByte: true);
    }

}
