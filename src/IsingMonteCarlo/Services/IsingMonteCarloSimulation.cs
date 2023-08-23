using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Services;

public sealed class IsingMonteCarloSimulation
{
    private const int LatticeSizeLowerBound = 3;

    private readonly List<List<int>> _neighboursIndices;
    private readonly IHamiltonian<int> _hamiltonian;
    private double _beta;
    private double _j;
    private double _h;
    private double? _jY;

    public IsingMonteCarloSimulation(
        int dimension,
        int latticeLength,
        IEnumerable<int>? initialSpinConfiguration = null)
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

        TotalSpinsCount = Convert.ToInt32(Math.Pow(latticeLength, dimension));

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

        Dimension = dimension;
        LatticeLength = latticeLength;
        Lattice = new NearestNeighbourNDIsingLattice<int>(dimension, latticeLength, initialSpins);

        _hamiltonian = new IsingHamiltonian(Lattice);
        _neighboursIndices = Lattice.NeighboursIndices;
        SpatialVectors = Lattice.SpatialVectors;
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public NearestNeighbourNDIsingLattice<int> Lattice { get; }

    public List<List<int>> SpatialVectors { get; }

    public void RunMonteCarlo(
        double beta,
        double j,
        double h,
        int? iterationLimit,
        SpinUpdateMethod spinUpdateMethod,
        double? jY = null,
        int? randomSeed = null)
    {
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
                _hamiltonian.FlipSpin(chosenSite);
            }
        }
    }

    private bool FlipWithMetropolis(double randomProbability, int siteIndex)
    {
        return randomProbability
               <= Math.Min(1, Math.Exp(-_beta * _hamiltonian.GetDeltaEnergyOfSite(siteIndex, _j, _h, _jY)));
    }

    private bool FlipWithGlauber(double randomProbability, int siteIndex)
    {
        return randomProbability
               <= 1 / (1 + Math.Exp(_beta * _hamiltonian.GetDeltaEnergyOfSite(siteIndex, _j, _h, _jY)));
    }

    private bool FlipSpinWithWolff(double randomProbability, int siteIndex)
    {
        throw new NotImplementedException();
    }
}
