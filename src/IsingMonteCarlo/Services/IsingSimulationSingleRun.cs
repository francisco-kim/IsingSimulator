using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;

namespace IsingMonteCarlo.Services;

public sealed class IsingSimulationSingleRun
{
    private const int ThermalisationStepsInLatticeSizeUnit = 20_000;

    private readonly SpinUpdateMethod _spinUpdateMethod;
    private readonly double _temperature;
    private readonly double _beta;
    private readonly double _j;
    private readonly double _h;
    private readonly double? _jY;
    private readonly int _dimension;
    private readonly int _latticeLength;
    private readonly int? _randomSeed;
    private readonly int _previousIterationCount;

    public IsingSimulationSingleRun(string? filename,
                                    int dimension,
                                    int latticeLength,
                                    double temperature,
                                    double j,
                                    double h,
                                    SpinUpdateMethod spinUpdateMethod = SpinUpdateMethod.Wolff,
                                    double? jY = null,
                                    double initialSpinDownRatio = 0.25,
                                    int? randomSeed = 41)
    {
        _temperature = temperature;
        _beta = 1 / temperature;
        _j = j;
        _h = h;
        _jY = jY;
        _dimension = dimension;
        _latticeLength = latticeLength;
        _randomSeed = randomSeed;
        _spinUpdateMethod = spinUpdateMethod;

        var (initialSpinConfiguration, _, previousIterationCount) =
            SpinConfigurationBuilder.InitialiseLattice(filename,
                                                       dimension,
                                                       latticeLength,
                                                       initialSpinDownRatio,
                                                       randomSeed);

        _previousIterationCount = previousIterationCount;

        Simulation = new IsingMonteCarloSimulation(dimension, latticeLength, initialSpinConfiguration);
    }

    public IsingMonteCarloSimulation Simulation { get; }

    public void RunWithMeasurements(int iterationStepsBetweenMeasurements,
                                    int measurementsCount,
                                    int thermalisationStepsInLatticeSizeUnit = ThermalisationStepsInLatticeSizeUnit,
                                    bool saveLattice = true)
    {
        Thermalise(thermalisationStepsInLatticeSizeUnit,
                   _spinUpdateMethod,
                   saveLattice);

        MeasurementsRun(iterationStepsBetweenMeasurements, measurementsCount);

        // if (saveLattice)
        // {
        //     LatticeConfigurationSaver.SaveLattice(Simulation.Lattice.Spins,
        //                                           _temperature,
        //                                           _previousIterationCount
        //                                           + thermalisationStepsInLatticeSizeUnit
        //                                           * Simulation.LatticeLength
        //                                           + measurementsCount
        //                                           * iterationStepsBetweenMeasurements,
        //                                           isBinary: false);
        // }
    }

    public void MeasurementsRun(int iterationStepsBetweenMeasurements, int measurementsCount)
    {
        Simulation.RunMonteCarloWithObservablesComputation(
            _beta,
            _j,
            _h,
            iterationStepsBetweenMeasurements,
            measurementsCount,
            _spinUpdateMethod,
            _randomSeed);

        Console.WriteLine($" M  = {Simulation.Magnetisation} +- {Simulation.MagnetisationSigma}");
        Console.WriteLine($"M^2 = {Simulation.MagnetisationSquared} +- {Simulation.MagnetisationSquaredSigma}");
        Console.WriteLine($"|M| = {Simulation.MagnetisationAbsolute} +- {Simulation.MagnetisationAbsoluteSigma}");
        Console.WriteLine($"Chi = {Simulation.Susceptibility}");
        Console.WriteLine($" Xi = {Simulation.RenormalisedCorrelationLength}\n");

        var susceptibilityList = Simulation.SusceptibilityList.Select(chi => $"{chi}");
        var renormalisedCorrelationLengthList = Simulation.RenormalisedCorrelationLengthList.Select(xi => $"{xi}");
        Console.WriteLine($"Chis = " + "{" + string.Join(", ", susceptibilityList) + "}");
        Console.WriteLine($" Xis = " + "{" + string.Join(", ", renormalisedCorrelationLengthList) + "}");
        // var correlationLengthListForPrintingX = correlationLengthList.Where(xi => xi.InXDirection is not double.NaN)
        //                                                              .Select(xi => $"{xi.InXDirection}");
        // var correlationLengthListForPrintingY = correlationLengthList.Where(xi => xi.InYDirection is not double.NaN)
        //                                                              .Select(xi => $"{xi.InYDirection}");
        // Console.WriteLine(
        //     string.Join(",", correlationLengthListForPrintingX.Concat(correlationLengthListForPrintingY).ToList()));
    }

    public void Thermalise(int thermalisationStepsInLatticeSizeUnit = ThermalisationStepsInLatticeSizeUnit,
                           SpinUpdateMethod spinUpdateMethod = SpinUpdateMethod.Wolff,
                           bool saveLattice = true)
    {
        Simulation.RunMonteCarlo(
            _beta,
            _j,
            _h,
            thermalisationStepsInLatticeSizeUnit * Simulation.TotalSpinsCount,
            spinUpdateMethod,
            _randomSeed);

        if (saveLattice)
        {
            LatticeConfigurationSaver.SaveLattice(Simulation.Lattice.Spins,
                                                  _temperature,
                                                  _previousIterationCount
                                                  + thermalisationStepsInLatticeSizeUnit
                                                  * Simulation.TotalSpinsCount,
                                                  isBinary: false);
        }

        Console.WriteLine($"M (T={_temperature}) = {Simulation.Hamiltonian.GetAverageMagnetisation(_j, _h)}"
            + " after thermalisation.\n");
    }
}