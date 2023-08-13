using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Services;

public sealed class MonteCarloSimulation
{
    private const int LatticeSizeLowerBound = 3;

    private readonly NearestNeighbourNDIsingLattice _lattice;
    private readonly List<List<int>> _neighboursIndices;

    private double? _jY;
    private double _j;
    private double _h;
    private double _beta;
    private List<int> Spins;

    public MonteCarloSimulation(
        int dimension,
        int latticeLength)
    {
        if (dimension < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(dimension),
                $"The dimension must be greater than 1, but {dimension} was given.");
        }

        if (latticeLength < LatticeSizeLowerBound)
        {
            throw new ArgumentOutOfRangeException(
                nameof(latticeLength),
                $"The lattice length be greater than {LatticeSizeLowerBound}, but {latticeLength} was given.");
        }

        Dimension = dimension;
        LatticeLength = latticeLength;
        TotalSpinsCount = Convert.ToInt32(Math.Pow(latticeLength, dimension));

        _lattice = new NearestNeighbourNDIsingLattice(dimension, latticeLength);
        _neighboursIndices = _lattice.NeighboursIndices;
        SpatialVectors = _lattice.SpatialVectors;

        Spins = Enumerable.Repeat(1, TotalSpinsCount).ToList();

        // var dimension = 3;
        // var sites = new List<List<int>>(latticeSize * latticeSize) { Enumerable.Repeat(0, dimension).ToList() };

        // foreach (var d in Enumerable.Range(0, dimension))
        // {
        // 	var sitesToAdd = new List<List<int>>(latticeSize - 1);
        // 	foreach(var siteToFill in sites)
        // 	{
        // 		for (int l = 1; l < latticeSize; l++)
        // 		{
        // 			var newSite = siteToFill.ToList();
        // 			newSite[d] = l;
        // 			sitesToAdd.Add(newSite);
        // 		}
        // 	}
        // 	sites.AddRange(sitesToAdd);
        // }

        // sites.OrderBy(x => x.FirstOrDefault());
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public List<List<int>> SpatialVectors { get; }

    // public List<int> Spins { get; }

    public void RunMonteCarlo(
        double beta,
        double j,
        double h,
        int? iterationLimit,
        SpinUpdateMethod spinUpdateMethod,
        IEnumerable<int>? initialSpinConfiguration = null,
        double? jY = null,
        int? randomSeed = null)
    {
        var initialSpins = initialSpinConfiguration is not null ?
                    initialSpinConfiguration.ToList() :
                    Enumerable.Repeat(1, TotalSpinsCount).ToList();

        if (initialSpins.Count != TotalSpinsCount)
        {
            throw new ArgumentException(
                $"There must be (lattice length)^(dimension) spins, but {initialSpins.Count} spins were given.",
                nameof(initialSpinConfiguration));
        }

        if (initialSpins.Any(spin  => (spin != 1) && (spin != -1)))
        {
            throw new ArgumentException(
                $"The spins must have integer values +1 or -1.",
                nameof(initialSpinConfiguration));
        }

        if (jY is not null && Dimension != 2)
        {
            throw new ArgumentException(
                "The coupling constant $J_{Y}$ is valid only in the 2D Ising model.",
                nameof(jY));
        }

        _beta = beta;
        _j = j;
        _h = h;
        _jY = jY;
        Spins = initialSpins;

        var random = randomSeed is not null ? new Random((int)randomSeed) : new Random();
        // bool flip = false;
        // switch (spinUpdateMethod)
        // {
        //     case SpinUpdateMethod.Metropolis:
        //     case SpinUpdateMethod.Glauber:
        //     case SpinUpdateMethod.Wolff:
        // }

        if (iterationLimit is null)
        {
            while (true)
            {
                flipSpinWithChosenUpdateMethod();
            }
        }
        else
        {
            int iterationCount = 0;
            while (iterationCount < iterationLimit)
            {
                flipSpinWithChosenUpdateMethod();
                ++iterationCount;
            }
        }

        void flipSpinWithChosenUpdateMethod()
        {
            var chosenSite = random.Next(0, TotalSpinsCount);
            var flip = FlipWithGlauber(randomProbability: random.NextDouble(), chosenSite);
            // var flip = FlipWithMetropolis(randomProbability: random.NextDouble(), chosenSite);

            if (flip)
            {
                FlipSpin(chosenSite);
            }
        }
    }

    public double GetTotalEnergy()
    {
        return 0.5 * Enumerable.Range(0, TotalSpinsCount)
                               .Select(GetEnergyOfSite)
                               .Sum();
    }

    /// <summary>
    ///     Calculates the average energy per site.
    ///     Its absolute value is the hypercube dimension D
    ///     (with an isotropic J-coupling and without an external field).
    /// </summary>
    public double GetAverageEnergy() => GetTotalEnergy() / TotalSpinsCount;

    public double GetAverageMagnetisation() => (double)Spins.Sum() / TotalSpinsCount;

    public double GetDeltaEnergyOfSite(int spinIndex) => -2.0 * GetEnergyOfSite(spinIndex);


    internal double GetEnergyOfSite(
        int spinIndex)
    {
        var spinValue = Spins[spinIndex];

        // Possible in 2D only
        if (_jY is not null && Dimension is 2)
        {
            var xBonds = _j * _neighboursIndices[spinIndex].Take(2)
                                                          .Select(
                                                            neighbourIndex => 
                                                                spinValue
                                                                * Spins[neighbourIndex])
                                                          .Sum();
            var yBonds = _jY * _neighboursIndices[spinIndex].Skip(2).Take(2)
                                                           .Select(
                                                            neighbourIndex => 
                                                                spinValue
                                                                * Spins[neighbourIndex])
                                                            .Sum();
            
            return xBonds + (double)yBonds + _h * spinValue;
        }

        return _j * _neighboursIndices[spinIndex].Select(neighbourIndex => Spins[neighbourIndex]
                                                                          * spinValue)
                                                .Sum()
               + _h * spinValue;
    }

    internal double FlipSpin(int spinIndex) => Spins[spinIndex] = -Spins[spinIndex];

    private bool FlipWithMetropolis(double randomProbability, int siteIndex)
    {
        return randomProbability
               <= Math.Min(1, Math.Exp(-_beta * GetDeltaEnergyOfSite(siteIndex)));
    }

    private bool FlipWithGlauber(double randomProbability, int siteIndex)
    {
        return randomProbability
               <= 1 / (1 + Math.Exp(_beta * GetDeltaEnergyOfSite(siteIndex)));
    }

    private bool FlipClusterWithWolff(double randomProbability, int siteIndex) => throw new NotImplementedException();
}
