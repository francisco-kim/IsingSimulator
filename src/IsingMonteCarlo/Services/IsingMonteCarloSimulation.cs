using System.ComponentModel;

using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;
using IsingMonteCarlo.Representations.SpinDynamics;

namespace IsingMonteCarlo.Services;

public sealed class IsingMonteCarloSimulation
{
    private const int LatticeSizeLowerBound = 3;

    private readonly List<List<int>> _neighboursIndices;
    private readonly IHamiltonian<int> _hamiltonian;

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

        var initialSpins = initialSpinConfiguration is not null
                               ? initialSpinConfiguration.ToList()
                               : Enumerable.Repeat(element: 1, TotalSpinsCount).ToList();

        if (initialSpins.Count != TotalSpinsCount)
        {
            throw new ArgumentException(
                $"There must be (lattice length)^(dimension) spins, but {initialSpins.Count} spins were given.",
                nameof(initialSpinConfiguration));
        }

        if (initialSpins.Any(spin => spin != 1 && spin != -1))
        {
            throw new ArgumentException(
                message: "The spins must have integer values +1 or -1.",
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

    public int TotalEnergy { get; set; }

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
                message: "The coupling constant $J_{Y}$ is valid only in the 2D Ising model.",
                nameof(jY));
        }

        if (!Enum.IsDefined(typeof(SpinUpdateMethod), spinUpdateMethod))
        {
            throw new InvalidEnumArgumentException(
                nameof(spinUpdateMethod),
                (int)spinUpdateMethod,
                typeof(SpinUpdateMethod));
        }

        ISpinDynamics spinDynamics = spinUpdateMethod switch
        {
            SpinUpdateMethod.Metropolis => new MetropolisDynamics(_hamiltonian, randomSeed),
            SpinUpdateMethod.Glauber => new GlauberDynamics(_hamiltonian, randomSeed),
            SpinUpdateMethod.Wolff => new WolffClusterDynamics(_hamiltonian, randomSeed),
            _ => new GlauberDynamics(_hamiltonian, randomSeed)
        };

        if (iterationLimit is null)
        {
            while (true)
            {
                spinDynamics.FlipSpin(beta, j, h, jY);
            }
        }

        var iterationCount = 0;
        while (iterationCount < iterationLimit)
        {
            spinDynamics.FlipSpin(beta, j, h, jY);
            ++iterationCount;
        }
    }
}
