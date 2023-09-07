using System.Numerics;

using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;
using IsingMonteCarlo.Representations.SpinDynamics;

namespace IsingMonteCarlo.Services;

public sealed class IsingMonteCarloSimulation
{
    private const int LatticeSizeLowerBound = 3;

    private readonly double _q1;
    private readonly double _q2;

    private int _measurementsCount;
    private double _magnetisationSum;
    private double _magnetisationSquaredSum;
    private double _magnetisationAbsoluteSum;
    private double _energySum;
    private double _structureFactorQ1XContributionFromRealSpinQ;
    private double _structureFactorQ1XContributionFromImaginarySpinQ;
    private double _structureFactorQ2XContributionFromRealSpinQ;
    private double _structureFactorQ2XContributionFromImaginarySpinQ;
    private double _structureFactorQ1YContributionFromRealSpinQ;
    private double _structureFactorQ1YContributionFromImaginarySpinQ;
    private double _structureFactorQ2YContributionFromRealSpinQ;
    private double _structureFactorQ2YContributionFromImaginarySpinQ;
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
        SpatialVectors = Lattice.SpatialVectors;

        _spinUpdateMethod = null;
       _spinDynamics = new GlauberDynamics(Hamiltonian);
 
        _q1 = 2.0 * Math.PI / Lattice.LatticeLength;
        _q2 = 2.0 * _q1;

        Magnetisation = double.NaN;
        MagnetisationSquared = double.NaN;
        MagnetisationAbsolute = double.NaN;
        Energy = double.NaN;
        Susceptibility = double.NaN;
        CorrelationLengthX = double.NaN;
        CorrelationLengthY = double.NaN;
        RenormalisedCorrelationLength = double.NaN;

        _measurementsCount = 0;
        MagnetisationList = new List<double>();
        MagnetisationSquaredList = new List<double>();
        MagnetisationAbsoluteList = new List<double>();
        EnergyList = new List<double>();
        _magnetisationSum = 0.0;
        _magnetisationSquaredSum = 0.0;
        _magnetisationAbsoluteSum = 0.0;
        _energySum = 0.0;
        _structureFactorQ1XContributionFromRealSpinQ = 0.0;
        _structureFactorQ1XContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ2XContributionFromRealSpinQ = 0.0;
        _structureFactorQ2XContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ1YContributionFromRealSpinQ = 0.0;
        _structureFactorQ1YContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ2YContributionFromRealSpinQ = 0.0;
        _structureFactorQ2YContributionFromImaginarySpinQ = 0.0;
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public NearestNeighbourNDIsingLattice<int> Lattice { get; }

    public IHamiltonian<int> Hamiltonian { get; }

    public List<List<int>> SpatialVectors { get; }

    public double Magnetisation { get; private set; }

    public double MagnetisationSquared { get; private set; }

    public double MagnetisationAbsolute { get; private set; }

    public double Energy { get; private set; }

    public double Susceptibility { get; private set; }

    public double CorrelationLengthX { get; private set; }

    public double CorrelationLengthY { get; private set; }

    public double RenormalisedCorrelationLength { get; private set; }

    public List<double> MagnetisationList { get; private set; }

    public List<double> MagnetisationSquaredList { get; private set; }

    public List<double> MagnetisationAbsoluteList { get; private set; }

    public List<double> EnergyList { get;  private set; }

    public void RunMonteCarlo(
        double beta,
        double j,
        double h,
        long? iterationLimit,
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
            _spinUpdateMethod = spinUpdateMethod;
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

        if (_spinUpdateMethod is SpinUpdateMethod.Wolff)
        {
            _spinDynamics.EmptyQueue(beta, j, h, jY);
        }
    }

    public void RunMonteCarloWithObservablesComputation(
        double beta,
        double j,
        double h,
        int iterationsNeededForSingleChiXiMeasurement,
        int measurementsCountForChiXiExpectationValue,
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

        if (iterationsNeededForSingleChiXiMeasurement <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(iterationsNeededForSingleChiXiMeasurement),
                "The observables xi and chi cannot be computed in less than two steps - the expectation value/average is needed.");
        }

        ResetObservables();

        var totalIterationCount = iterationsNeededForSingleChiXiMeasurement * measurementsCountForChiXiExpectationValue;
        var iterationCount = 0;

        if (measurementsCountForChiXiExpectationValue is 0)
        {
            iterateWithObservablesMeasurements(isInfiniteLoop: true);
        }

        iterateWithObservablesMeasurements(isInfiniteLoop: false);

        void iterateWithObservablesMeasurements(bool isInfiniteLoop)
        {
            var iterationConditionIsSatisfied = true;
            while (iterationConditionIsSatisfied)
            {
                _spinDynamics.FlipSpin(beta, j, h, jY);

                if (iterationCount % iterationsNeededForSingleChiXiMeasurement == 0)
                {
                    if (_spinUpdateMethod is SpinUpdateMethod.Wolff)
                    {
                        _spinDynamics.EmptyQueue(beta, j, h, jY);
                    }

                    var magnetisationMeasurement = Hamiltonian.GetAverageMagnetisation(j, h, jY);
                    var magnetisationSquaredMeasurement = magnetisationMeasurement * magnetisationMeasurement;
                    var magnetisationAbsoluteMeasurement = Math.Abs(magnetisationMeasurement);
                    var energyMeasurement = Hamiltonian.GetAverageEnergy(j, h, jY);

                    MagnetisationList.Add(magnetisationMeasurement);
                    MagnetisationSquaredList.Add(magnetisationSquaredMeasurement);
                    MagnetisationAbsoluteList.Add(magnetisationAbsoluteMeasurement);
                    EnergyList.Add(energyMeasurement);

                    _magnetisationSum += magnetisationMeasurement;
                    _magnetisationSquaredSum += magnetisationSquaredMeasurement;
                    _magnetisationAbsoluteSum += magnetisationAbsoluteMeasurement;
                    _energySum += energyMeasurement;
                    MeasureStructureFactorSQ1AndSQ2();

                    _measurementsCount++;

                    Magnetisation = _magnetisationSum / _measurementsCount;
                    MagnetisationSquared = _magnetisationSquaredSum / _measurementsCount;
                    MagnetisationAbsolute = _magnetisationAbsoluteSum / _measurementsCount;
                    Energy = _energySum / _measurementsCount;

                    Susceptibility = TotalSpinsCount
                                   * beta
                                   * (MagnetisationSquared - MagnetisationAbsolute * MagnetisationAbsolute);
                    (CorrelationLengthX, CorrelationLengthY) = GetCorrelationLengthInXYDirections();

                    if (!double.IsNaN(CorrelationLengthX) && !double.IsNaN(CorrelationLengthX))
                    {
                        RenormalisedCorrelationLength = (CorrelationLengthX + CorrelationLengthY) / 2.0 / LatticeLength;
                    }
                    else if (!double.IsNaN(CorrelationLengthX))
                    {
                        RenormalisedCorrelationLength = CorrelationLengthX / LatticeLength;
                    }
                    else if (!double.IsNaN(CorrelationLengthY))
                    {
                        RenormalisedCorrelationLength = CorrelationLengthY / LatticeLength;
                    }
                }

                ++iterationCount;

                if (isInfiniteLoop)
                {
                    if (iterationCount >= int.MaxValue || measurementsCountForChiXiExpectationValue >= 1000)
                    {
                        iterationCount = 0;
                    }
                }
                else
                {
                    iterationConditionIsSatisfied = iterationCount < totalIterationCount;
                }
            }
        }
    }

    internal void ResetObservables()
    {
        Magnetisation = 0.0;
        MagnetisationSquared = 0.0;
        MagnetisationAbsolute = 0.0;
        Energy = 0.0;
        Susceptibility = 0.0;
        CorrelationLengthX = 0.0;
        CorrelationLengthY = 0.0;
        RenormalisedCorrelationLength = 0.0;

        _measurementsCount = 0;
        MagnetisationList = new List<double>();
        MagnetisationSquaredList = new List<double>();
        MagnetisationAbsoluteList = new List<double>();
        EnergyList = new List<double>();
        _magnetisationSum = 0.0;
        _magnetisationSquaredSum = 0.0;
        _magnetisationAbsoluteSum = 0.0;
        _energySum = 0.0;
        _structureFactorQ1XContributionFromRealSpinQ = 0.0;
        _structureFactorQ1XContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ2XContributionFromRealSpinQ = 0.0;
        _structureFactorQ2XContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ1YContributionFromRealSpinQ = 0.0;
        _structureFactorQ1YContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ2YContributionFromRealSpinQ = 0.0;
        _structureFactorQ2YContributionFromImaginarySpinQ = 0.0;
    }

    // Only if coupling strength is symmetrical in all directions
    private (double InXDirection, double InYDirection) GetCorrelationLengthInXYDirections()
    {
        var sQ1X = (_structureFactorQ1XContributionFromRealSpinQ + _structureFactorQ1XContributionFromImaginarySpinQ)
                / _measurementsCount;
        var sQ2X = (_structureFactorQ2XContributionFromRealSpinQ + _structureFactorQ2XContributionFromImaginarySpinQ)
                / _measurementsCount;

        var sQ1Y = (_structureFactorQ1YContributionFromRealSpinQ + _structureFactorQ1YContributionFromImaginarySpinQ)
                / _measurementsCount;
        var sQ2Y = (_structureFactorQ2YContributionFromRealSpinQ + _structureFactorQ2YContributionFromImaginarySpinQ)
                / _measurementsCount;

        var alternativeCorrelationLengthX = 1.0 / _q1 * Math.Sqrt((sQ1X / sQ2X - 1.0) / (4.0 - sQ1X / sQ2X));
        var alternativeCorrelationLengthY = 1.0 / _q1 * Math.Sqrt((sQ1Y / sQ2Y - 1.0) / (4.0 - sQ1Y / sQ2Y));

        var factor = 1.0 / Math.Sqrt((1.0 + Dimension) * (3.0 + Dimension) / (8.0 * Dimension));
        return (alternativeCorrelationLengthX * factor, alternativeCorrelationLengthY * factor);
    }

    private void MeasureStructureFactorSQ1AndSQ2()
    {
        var (realQ1X, imaginaryQ1X) = GetFTSpinQInOneDirection(_q1, 0);
        var (realQ2X, imaginaryQ2X) = GetFTSpinQInOneDirection(_q2, 0);

        var (realQ1Y, imaginaryQ1Y) = GetFTSpinQInOneDirection(_q1, 1);
        var (realQ2Y, imaginaryQ2Y) = GetFTSpinQInOneDirection(_q2, 1);

        // In the x-direction
        _structureFactorQ1XContributionFromRealSpinQ += realQ1X * realQ1X;
        _structureFactorQ1XContributionFromImaginarySpinQ += imaginaryQ1X * imaginaryQ1X;
        _structureFactorQ2XContributionFromRealSpinQ += realQ2X * realQ2X;
        _structureFactorQ2XContributionFromImaginarySpinQ += imaginaryQ2X * imaginaryQ2X;

        // In the y-direction
        _structureFactorQ1YContributionFromRealSpinQ += realQ1Y * realQ1Y;
        _structureFactorQ1YContributionFromImaginarySpinQ += imaginaryQ1Y * imaginaryQ1Y;
        _structureFactorQ2YContributionFromRealSpinQ += realQ2Y * realQ2Y;
        _structureFactorQ2YContributionFromImaginarySpinQ += imaginaryQ2Y * imaginaryQ2Y;
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
