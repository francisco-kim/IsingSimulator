using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations.SpinDynamics;

namespace IsingMonteCarlo.Representations;

public sealed class IsingMonteCarloSimulator
{
    private const int LatticeSizeLowerBound = 3;

    private readonly double _q1;
    private readonly double _q2;
    private readonly double[] _cosQ1;
    private readonly double[] _sinQ1;
    private readonly double[] _cosQ2;
    private readonly double[] _sinQ2;

    private int _measurementsCount;
    private double _magnetisationSum;
    private double _magnetisationSquaredSum;
    private double _magnetisationAbsoluteSum;
    private double _energySum;
    private double _energySquaredSum;
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

    public IsingMonteCarloSimulator(
        int dimension,
        int latticeLength,
        IEnumerable<int>? initialSpinConfiguration = null)
    {
        if (dimension < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(dimension),
                $"The dimension must be at least 1, but {dimension} was given.");
        }

        if (latticeLength < LatticeSizeLowerBound)
        {
            throw new ArgumentOutOfRangeException(
                nameof(latticeLength),
                $"The lattice length must be at least {LatticeSizeLowerBound}, but {latticeLength} was given.");
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

        _cosQ1 = new double[latticeLength];
        _sinQ1 = new double[latticeLength];
        _cosQ2 = new double[latticeLength];
        _sinQ2 = new double[latticeLength];
        for (var position = 0; position < latticeLength; position++)
        {
            _cosQ1[position] = Math.Cos(_q1 * position);
            _sinQ1[position] = Math.Sin(_q1 * position);
            _cosQ2[position] = Math.Cos(_q2 * position);
            _sinQ2[position] = Math.Sin(_q2 * position);
        }

        ResetObservables();
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public NearestNeighbourNDIsingLattice<int> Lattice { get; }

    public IHamiltonian<int> Hamiltonian { get; }

    public List<List<int>> SpatialVectors { get; }

    /// <summary>
    ///     The sites of the Wolff cluster currently being grown, or null when
    ///     the active dynamics is not Wolff. Intended for visualisation.
    /// </summary>
    public IReadOnlyList<int>? CurrentWolffCluster =>
        (_spinDynamics as WolffClusterDynamics)?.CurrentClusterSites;

    /// <summary>
    ///     True while a Wolff cluster is partially grown.
    /// </summary>
    public bool HasPendingWolffCluster =>
        _spinDynamics is WolffClusterDynamics { HasPendingCluster: true };

    /// <summary>
    ///     Clears all measurement accumulators so that a new measurement run
    ///     is statistically independent of the previous ones.
    /// </summary>
    public void ResetObservables()
    {
        _measurementsCount = 0;
        _magnetisationSum = 0.0;
        _magnetisationSquaredSum = 0.0;
        _magnetisationAbsoluteSum = 0.0;
        _energySum = 0.0;
        _energySquaredSum = 0.0;
        _structureFactorQ1XContributionFromRealSpinQ = 0.0;
        _structureFactorQ1XContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ2XContributionFromRealSpinQ = 0.0;
        _structureFactorQ2XContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ1YContributionFromRealSpinQ = 0.0;
        _structureFactorQ1YContributionFromImaginarySpinQ = 0.0;
        _structureFactorQ2YContributionFromRealSpinQ = 0.0;
        _structureFactorQ2YContributionFromImaginarySpinQ = 0.0;
    }

    public void RunMonteCarlo(
        double beta,
        double j,
        double h,
        long? iterationLimit,
        SpinUpdateMethod spinUpdateMethod,
        double? jY = null,
        int? randomSeed = null)
    {
        ValidateParameters(spinUpdateMethod, h, jY);
        EnsureSpinDynamics(spinUpdateMethod, randomSeed);

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
            _spinDynamics.EmptyQueue(beta, j, h, jY, verbose: true);
        }
    }

    /// <summary>
    ///     Runs a fixed number of spin-update steps without flushing a pending
    ///     Wolff cluster afterwards, so callers (e.g. an interactive renderer)
    ///     may observe the lattice mid-cluster-growth. Call
    ///     <see cref="FlushWolffQueue" /> before taking measurements or
    ///     switching to another update method.
    /// </summary>
    public void RunSteps(
        long stepCount,
        double beta,
        double j,
        double h,
        SpinUpdateMethod spinUpdateMethod,
        double? jY = null,
        int? randomSeed = null)
    {
        ValidateParameters(spinUpdateMethod, h, jY);
        EnsureSpinDynamics(spinUpdateMethod, randomSeed);

        for (long step = 0; step < stepCount; step++)
        {
            _spinDynamics.FlipSpin(beta, j, h, jY);
        }
    }

    /// <summary>
    ///     Completes any partially grown Wolff cluster. No-op for the local
    ///     update methods.
    /// </summary>
    public void FlushWolffQueue(
        double beta,
        double j,
        double h,
        double? jY = null)
    {
        if (_spinUpdateMethod is SpinUpdateMethod.Wolff)
        {
            _spinDynamics.EmptyQueue(beta, j, h, jY, verbose: false);
        }
    }

    public BasicObservables RunMonteCarloWithObservablesComputation(
        double beta,
        double j,
        double h,
        int iterationsNeededForSingleChiXiMeasurement,
        int measurementsCountForChiXiExpectationValue,
        SpinUpdateMethod spinUpdateMethod,
        double? jY = null,
        int? randomSeed = null)
    {
        ValidateParameters(spinUpdateMethod, h, jY);
        EnsureSpinDynamics(spinUpdateMethod, randomSeed);

        if (iterationsNeededForSingleChiXiMeasurement <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(iterationsNeededForSingleChiXiMeasurement),
                "The observables xi and chi cannot be computed in less than two steps - the expectation value/average is needed.");
        }

        double magnetisation = double.NaN;
        double magnetisationSquared = double.NaN;
        double magnetisationAbsolute = double.NaN;
        double energy = double.NaN;
        double susceptibility = double.NaN;
        double specificHeat = double.NaN;
        double correlationLengthX = double.NaN;
        double correlationLengthY = double.NaN;
        double renormalisedCorrelationLength = double.NaN;
        List<double> magnetisationList = new();
        List<double> magnetisationSquaredList = new();
        List<double> magnetisationAbsoluteList = new();
        List<double> energyList = new();

        var totalIterationCount =
            (long)iterationsNeededForSingleChiXiMeasurement * measurementsCountForChiXiExpectationValue;
        var iterationCount = 0L;
        while (iterationCount < totalIterationCount)
        {
            _spinDynamics.FlipSpin(beta, j, h, jY);

            // Sample at the end of each decorrelation block, not at its start.
            if ((iterationCount + 1) % iterationsNeededForSingleChiXiMeasurement == 0)
            {
                if (_spinUpdateMethod is SpinUpdateMethod.Wolff)
                {
                    _spinDynamics.EmptyQueue(beta, j, h, jY, verbose: true);
                }

                var magnetisationMeasurement = Hamiltonian.GetAverageMagnetisation(j, h, jY);
                var magnetisationSquaredMeasurement = magnetisationMeasurement * magnetisationMeasurement;
                var magnetisationAbsoluteMeasurement = Math.Abs(magnetisationMeasurement);
                var energyMeasurement = Hamiltonian.GetAverageEnergy(j, h, jY);

                magnetisationList.Add(magnetisationMeasurement);
                magnetisationSquaredList.Add(magnetisationSquaredMeasurement);
                magnetisationAbsoluteList.Add(magnetisationAbsoluteMeasurement);
                energyList.Add(energyMeasurement);

                _magnetisationSum += magnetisationMeasurement;
                _magnetisationSquaredSum += magnetisationSquaredMeasurement;
                _magnetisationAbsoluteSum += magnetisationAbsoluteMeasurement;
                _energySum += energyMeasurement;
                _energySquaredSum += energyMeasurement * energyMeasurement;
                MeasureStructureFactorSQ1AndSQ2();

                _measurementsCount++;

                magnetisation = _magnetisationSum / _measurementsCount;
                magnetisationSquared = _magnetisationSquaredSum / _measurementsCount;
                magnetisationAbsolute = _magnetisationAbsoluteSum / _measurementsCount;
                energy = _energySum / _measurementsCount;
                var energySquared = _energySquaredSum / _measurementsCount;

                susceptibility = TotalSpinsCount
                               * beta
                               * (magnetisationSquared - magnetisationAbsolute * magnetisationAbsolute);
                specificHeat = TotalSpinsCount
                             * beta
                             * beta
                             * (energySquared - energy * energy);
                (correlationLengthX, correlationLengthY) = GetCorrelationLengthInXYDirections();

                if (!double.IsNaN(correlationLengthX) && !double.IsNaN(correlationLengthY))
                {
                    renormalisedCorrelationLength = (correlationLengthX + correlationLengthY) / 2.0 / LatticeLength;
                }
                else if (!double.IsNaN(correlationLengthX))
                {
                    renormalisedCorrelationLength = correlationLengthX / LatticeLength;
                }
                else if (!double.IsNaN(correlationLengthY))
                {
                    renormalisedCorrelationLength = correlationLengthY / LatticeLength;
                }
            }

            ++iterationCount;
        }

        return new BasicObservables(
            magnetisation,
            magnetisationSquared,
            magnetisationAbsolute,
            energy,
            correlationLengthX,
            correlationLengthY,
            renormalisedCorrelationLength,
            susceptibility,
            specificHeat,
            magnetisationList,
            magnetisationSquaredList,
            magnetisationAbsoluteList,
            energyList);
    }

    private void ValidateParameters(
        SpinUpdateMethod spinUpdateMethod,
        double h,
        double? jY)
    {
        if (jY is not null && Dimension != 2)
        {
            throw new ArgumentException(
                "The coupling constant $J_{Y}$ is valid only in the 2D Ising model.",
                nameof(jY));
        }

        if (h is not 0.0 && spinUpdateMethod is SpinUpdateMethod.Wolff)
        {
            throw new ArgumentException(
                "The Wolff single-cluster algorithm is not allowed "
              + $"if the external field is present ({nameof(h)} = {h} here).",
                nameof(h));
        }
    }

    private void EnsureSpinDynamics(SpinUpdateMethod spinUpdateMethod, int? randomSeed)
    {
        if (spinUpdateMethod == _spinUpdateMethod)
        {
            return;
        }

        // Note: the seed only applies when the dynamics object is (re)created,
        // i.e. on the first run or when the update method changes.
        _spinDynamics = spinUpdateMethod switch
        {
            SpinUpdateMethod.Metropolis => new MetropolisDynamics(Hamiltonian, randomSeed),
            SpinUpdateMethod.Glauber => new GlauberDynamics(Hamiltonian, randomSeed),
            SpinUpdateMethod.Wolff => new WolffClusterDynamics(Hamiltonian, randomSeed),
            _ => new GlauberDynamics(Hamiltonian, randomSeed)
        };
        _spinUpdateMethod = spinUpdateMethod;
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

        // Two-momentum estimator from the Ornstein-Zernike form S(q) ~ 1/(q^2 + xi^-2)
        // evaluated at q1 = 2*pi/L and q2 = 2*q1.
        var correlationLengthX = 1.0 / _q1 * Math.Sqrt((sQ1X / sQ2X - 1.0) / (4.0 - sQ1X / sQ2X));
        var correlationLengthY = 1.0 / _q1 * Math.Sqrt((sQ1Y / sQ2Y - 1.0) / (4.0 - sQ1Y / sQ2Y));

        return (correlationLengthX, correlationLengthY);
    }

    private void MeasureStructureFactorSQ1AndSQ2()
    {
        var (realQ1X, imaginaryQ1X) = GetFTSpinQInOneDirection(useQ2: false, componentIndex: 0);
        var (realQ2X, imaginaryQ2X) = GetFTSpinQInOneDirection(useQ2: true, componentIndex: 0);

        // In the x-direction
        _structureFactorQ1XContributionFromRealSpinQ += realQ1X * realQ1X;
        _structureFactorQ1XContributionFromImaginarySpinQ += imaginaryQ1X * imaginaryQ1X;
        _structureFactorQ2XContributionFromRealSpinQ += realQ2X * realQ2X;
        _structureFactorQ2XContributionFromImaginarySpinQ += imaginaryQ2X * imaginaryQ2X;

        if (Dimension < 2)
        {
            return;
        }

        var (realQ1Y, imaginaryQ1Y) = GetFTSpinQInOneDirection(useQ2: false, componentIndex: 1);
        var (realQ2Y, imaginaryQ2Y) = GetFTSpinQInOneDirection(useQ2: true, componentIndex: 1);

        // In the y-direction
        _structureFactorQ1YContributionFromRealSpinQ += realQ1Y * realQ1Y;
        _structureFactorQ1YContributionFromImaginarySpinQ += imaginaryQ1Y * imaginaryQ1Y;
        _structureFactorQ2YContributionFromRealSpinQ += realQ2Y * realQ2Y;
        _structureFactorQ2YContributionFromImaginarySpinQ += imaginaryQ2Y * imaginaryQ2Y;
    }

    private (double Real, double Imaginary) GetFTSpinQInOneDirection(bool useQ2, int componentIndex)
    {
        var cosTable = useQ2 ? _cosQ2 : _cosQ1;
        var sinTable = useQ2 ? _sinQ2 : _sinQ1;
        var spins = Lattice.Spins;

        var real = 0.0;
        var imaginary = 0.0;
        for (var spinIndex = 0; spinIndex < TotalSpinsCount; spinIndex++)
        {
            var position = SpatialVectors[spinIndex][componentIndex];

            // spin * exp(-i * q * position)
            real += spins[spinIndex] * cosTable[position];
            imaginary -= spins[spinIndex] * sinTable[position];
        }

        var factor = 1.0 / Math.Sqrt(TotalSpinsCount);

        return (factor * real, factor * imaginary);
    }
}
