using System.Numerics;

using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;
using IsingMonteCarlo.Representations.SpinDynamics;

namespace IsingMonteCarlo.Services;

public sealed class IsingMonteCarloSimulation
{
    private const int LatticeSizeLowerBound = 3;

    private readonly List<List<int>> _neighboursIndices;
    private readonly double _q1;
    private readonly double _q2;

    private SpinUpdateMethod? _spinUpdateMethod;
    private ISpinDynamics _spinDynamics;

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

        Hamiltonian = new IsingHamiltonian(Lattice);
        _neighboursIndices = Lattice.NeighboursIndices;
        SpatialVectors = Lattice.SpatialVectors;

        _spinUpdateMethod = null;
       _spinDynamics = new GlauberDynamics(Hamiltonian);
 
        _q1 = 2.0 * Math.PI / Lattice.LatticeLength;
        _q2 = 2.0 * _q1;

        StructureFactorQ1XContributionFromRealSpinQ = new List<double>();
        StructureFactorQ1XContributionFromImaginarySpinQ = new List<double>();
        StructureFactorQ2XContributionFromRealSpinQ = new List<double>();
        StructureFactorQ2XContributionFromImaginarySpinQ = new List<double>();
        StructureFactorQ1YContributionFromRealSpinQ = new List<double>();
        StructureFactorQ1YContributionFromImaginarySpinQ = new List<double>();
        StructureFactorQ2YContributionFromRealSpinQ = new List<double>();
        StructureFactorQ2YContributionFromImaginarySpinQ = new List<double>();
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public NearestNeighbourNDIsingLattice<int> Lattice { get; }

    public IHamiltonian<int> Hamiltonian { get; }

    public List<List<int>> SpatialVectors { get; }

    public List<double> StructureFactorQ1XContributionFromRealSpinQ { get; private set; }

    public List<double> StructureFactorQ1XContributionFromImaginarySpinQ { get; private set; }

    public List<double> StructureFactorQ2XContributionFromRealSpinQ { get; private set; }

    public List<double> StructureFactorQ2XContributionFromImaginarySpinQ { get; private set; }

    public List<double> StructureFactorQ1YContributionFromRealSpinQ { get; private set; }

    public List<double> StructureFactorQ1YContributionFromImaginarySpinQ { get; private set; }

    public List<double> StructureFactorQ2YContributionFromRealSpinQ { get; private set; }

    public List<double> StructureFactorQ2YContributionFromImaginarySpinQ { get; private set; }

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

        if (h is not 0.0 && spinUpdateMethod is SpinUpdateMethod.Wolff)
        {
            throw new ArgumentException("The Wolff single-cluster algorithm is not allowed "
                + $"if the extrnal field is present ({nameof(h)} = {h} here).", nameof(h));
        }

        if (spinUpdateMethod != _spinUpdateMethod)
        {
            _spinDynamics = spinUpdateMethod switch
            {
                SpinUpdateMethod.Metropolis => new MetropolisDynamics(Hamiltonian, randomSeed),
                SpinUpdateMethod.Glauber => new GlauberDynamics(Hamiltonian, randomSeed),
                SpinUpdateMethod.Wolff => new WolffClusterDynamics(Hamiltonian, randomSeed),
                _ => new GlauberDynamics(Hamiltonian, randomSeed)
            };
        }

        if (iterationLimit is null)
        {
            while (true)
            {
                _spinDynamics.FlipSpin(beta, j, h, jY);
            }
        }

        var iterationCount = 0;
        while (iterationCount < iterationLimit)
        {
            _spinDynamics.FlipSpin(beta, j, h, jY);
            ++iterationCount;
        }
    }

    // Only if coupling strength is symmetrical in all directions
    public (double InXDirection, double InYDirection) GetCorrelationLength()
    {
        MeasureStructureFactorSQ1AndSQ2();

        var sQ1X = (StructureFactorQ1XContributionFromRealSpinQ.Sum() + StructureFactorQ1XContributionFromImaginarySpinQ.Sum())
                / StructureFactorQ1XContributionFromRealSpinQ.Count;
        var sQ2X = (StructureFactorQ2XContributionFromRealSpinQ.Sum() + StructureFactorQ2XContributionFromImaginarySpinQ.Sum())
                / StructureFactorQ2XContributionFromRealSpinQ.Count;

        var sQ1Y = (StructureFactorQ1YContributionFromRealSpinQ.Sum() + StructureFactorQ1YContributionFromImaginarySpinQ.Sum())
                / StructureFactorQ1YContributionFromRealSpinQ.Count;
        var sQ2Y = (StructureFactorQ2YContributionFromRealSpinQ.Sum() + StructureFactorQ2YContributionFromImaginarySpinQ.Sum())
                / StructureFactorQ2YContributionFromRealSpinQ.Count;

        var alternativeCorrelationLengthX = 1.0 / _q1 * Math.Sqrt((sQ1X / sQ2X - 1.0) / (4.0 - sQ1X / sQ2X));
        var alternativeCorrelationLengthY = 1.0 / _q1 * Math.Sqrt((sQ1Y / sQ2Y - 1.0) / (4.0 - sQ1Y / sQ2Y));

        //return alternativeCorrelationLength
        //     / Math.Sqrt((1.0 + Dimension) * (3.0 + Dimension) / (8.0 * Dimension));
        return (alternativeCorrelationLengthX / LatticeLength, alternativeCorrelationLengthY / LatticeLength);
    }

    private void MeasureStructureFactorSQ1AndSQ2()
    {
        var (realQ1X, imaginaryQ1X) = GetFTSpinQInOneDirection(_q1, 0);
        var (realQ2X, imaginaryQ2X) = GetFTSpinQInOneDirection(_q2, 0);

        var (realQ1Y, imaginaryQ1Y) = GetFTSpinQInOneDirection(_q1, 1);
        var (realQ2Y, imaginaryQ2Y) = GetFTSpinQInOneDirection(_q2, 1);

        // In the x-direction
        StructureFactorQ1XContributionFromRealSpinQ.Add(realQ1X * realQ1X);
        StructureFactorQ1XContributionFromImaginarySpinQ.Add(imaginaryQ1X * imaginaryQ1X);
        StructureFactorQ2XContributionFromRealSpinQ.Add(realQ2X * realQ2X);
        StructureFactorQ2XContributionFromImaginarySpinQ.Add(imaginaryQ2X * imaginaryQ2X);

        // In the y-direction
        StructureFactorQ1YContributionFromRealSpinQ.Add(realQ1Y * realQ1Y);
        StructureFactorQ1YContributionFromImaginarySpinQ.Add(imaginaryQ1Y * imaginaryQ1Y);
        StructureFactorQ2YContributionFromRealSpinQ.Add(realQ2Y * realQ2Y);
        StructureFactorQ2YContributionFromImaginarySpinQ.Add(imaginaryQ2Y * imaginaryQ2Y);
    }

    private (double Real, double Imaginary) GetFTSpinQInOneDirection(double qX, int componentIndex)
    {
        var summationElements = Lattice.Spins.Select(
                                           (s, i) => s * Complex.Exp(-new Complex(0.0, 1.0) * qX * Lattice.SpatialVectors[i][componentIndex]))
                                       .ToList();

        var factor = 1.0 / Math.Sqrt(TotalSpinsCount);

        var real = factor * summationElements.Select(element => element.Real).Sum();
        var imaginary = factor * summationElements.Select(element => element.Imaginary).Sum();

        return (real, imaginary);
    }

    private (double Real, double Imaginary) GetFTSpinQ(List<double> q)
    {
        //if (q.Count != _dimension)
        //{
        //    throw new ArgumentException(
        //        $"The q-vector must be of the same dimension as the lattice, but it has dimension {q.Count}.",
        //        nameof(q));
        //}

        var summationElements = Lattice.Spins.Select(
                                           (s, i) => s
                                                   * Complex.Exp(
                                                         -new Complex(0.0, 1.0)
                                                       * new Complex(
                                                             scalarProduct(q, Lattice.SpatialVectors[i]),
                                                             0.0)))
                                       .ToList();

        var factor = 1.0 / Math.Sqrt(TotalSpinsCount);

        var real = factor * summationElements.Select(element => element.Real).Sum();
        var imaginary = factor * summationElements.Select(element => element.Real).Sum();

        return (real, imaginary);

        static double scalarProduct(List<double> v1, List<int> v2) =>
            v1.Zip(v2, (firstV, secondV) => firstV * secondV).Sum();
    }
}
