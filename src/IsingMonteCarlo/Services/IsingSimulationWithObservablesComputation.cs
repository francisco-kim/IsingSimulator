using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Services;

public sealed class IsingSimulationWithObservablesComputation
{
    private const int ThermalisationStepsInMCSweepUnit = 20_000;

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

    public IsingSimulationWithObservablesComputation(
        string? filename,
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
            SpinConfigurationBuilder.InitialiseLattice(
                filename,
                dimension,
                latticeLength,
                initialSpinDownRatio,
                randomSeed);

        _previousIterationCount = previousIterationCount;

        Simulation = new IsingMonteCarloSimulator(dimension, latticeLength, initialSpinConfiguration);
    }

    public IsingMonteCarloSimulator Simulation { get; }

    public Observables RunWithMeasurements(
        int iterationsNeededForSingleChiXiMeasurement = 5000,
        int measurementsCountForChiXiExpectationValue = 200,
        int measurementsRepetitionCountForChiXiVariance = 5,
        int thermalisationStepsInMCSweepUnit = ThermalisationStepsInMCSweepUnit,
        bool saveLattice = true,
        bool saveMeasurements = true,
        bool resetIterationCountDuringSave = false)
    {
        Thermalise(
            thermalisationStepsInMCSweepUnit,
            _spinUpdateMethod,
            saveLattice,
            resetIterationCountDuringSave);

        return MeasurementsRun(
            iterationsNeededForSingleChiXiMeasurement,
            measurementsCountForChiXiExpectationValue,
            measurementsRepetitionCountForChiXiVariance,
            saveMeasurements);

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

    public Observables MeasurementsRun(
        int iterationNeededForSingleChiXiMeasurement,
        int measurementsCountForChiXiExpectationValue,
        int measurementsRepetitionCountForChiXiVariance,
        bool saveMeasurements = true)
    {
        List<double> magnetisationList = new();
        List<double> magnetisationSquaredList = new();
        List<double> magnetisationAbsoluteList = new();
        List<double> energyList = new();
        List<double> susceptibilityList = new();
        List<double> renormalisedCorrelationLengthList = new();
        List<double> correlationLengthXList = new();
        List<double> correlationLengthYList = new();

        for (var measurementRepetition = 0;
             measurementRepetition < measurementsRepetitionCountForChiXiVariance;
             measurementRepetition++)
        {
            var basicObservables = Simulation.RunMonteCarloWithObservablesComputation(
                _beta,
                _j,
                _h,
                iterationNeededForSingleChiXiMeasurement,
                measurementsCountForChiXiExpectationValue,
                _spinUpdateMethod,
                _jY,
                _randomSeed);

            magnetisationList.AddRange(basicObservables.MagnetisationList);
            magnetisationSquaredList.AddRange(basicObservables.MagnetisationSquaredList);
            magnetisationAbsoluteList.AddRange(basicObservables.MagnetisationAbsoluteList);
            energyList.AddRange(basicObservables.EnergyList);

            susceptibilityList.Add(basicObservables.Susceptibility);

            if (!double.IsNaN(basicObservables.CorrelationLengthX))
            {
                correlationLengthXList.Add(basicObservables.CorrelationLengthX);
            }

            if (!double.IsNaN(basicObservables.CorrelationLengthY))
            {
                correlationLengthYList.Add(basicObservables.CorrelationLengthY);
            }

            if (!double.IsNaN(basicObservables.RenormalisedCorrelationLength))
            {
                renormalisedCorrelationLengthList.Add(basicObservables.RenormalisedCorrelationLength);
            }
        }

        var magnetisation = magnetisationList.Sum() / magnetisationList.Count;
        var magnetisationSigma = Math.Sqrt(GetVariance(magnetisationList, magnetisation));

        var magnetisationSquared = magnetisationSquaredList.Sum() / magnetisationSquaredList.Count;
        var magnetisationSquaredSigma = Math.Sqrt(GetVariance(magnetisationSquaredList, magnetisationSquared));

        var magnetisationAbsolute = magnetisationAbsoluteList.Sum() / magnetisationAbsoluteList.Count;
        var magnetisationAbsoluteSigma = Math.Sqrt(GetVariance(magnetisationAbsoluteList, magnetisationAbsolute));

        var energy = energyList.Sum() / energyList.Count;
        var energySigma = Math.Sqrt(GetVariance(energyList, energy));

        var correlationLengthX =
            correlationLengthXList.Sum() / correlationLengthXList.Count;
        var correlationLengthXSigma =
            Math.Sqrt(GetVariance(correlationLengthXList, correlationLengthX));
        var correlationLengthY =
            correlationLengthYList.Sum() / correlationLengthYList.Count;
        var correlationLengthYSigma =
            Math.Sqrt(GetVariance(correlationLengthYList, correlationLengthY));
        var renormalisedCorrelationLength =
            renormalisedCorrelationLengthList.Sum() / renormalisedCorrelationLengthList.Count;
        var renormalisedCorrelationLengthSigma =
            Math.Sqrt(GetVariance(renormalisedCorrelationLengthList, renormalisedCorrelationLength));

        var susceptibility = susceptibilityList.Sum() / susceptibilityList.Count;
        var susceptibilitySigma = Math.Sqrt(GetVariance(susceptibilityList, susceptibility));

        var magnetisationString = $" M  = {magnetisation} +- {magnetisationSigma}";
        var magnetisationSquaredString = $"M^2 = {magnetisationSquared} +- {magnetisationSquaredSigma}";
        var magnetisationAbsoluteString = $"|M| = {magnetisationAbsolute} +- {magnetisationAbsoluteSigma}";
        var energyString = $"E  = {energy} +- {energySigma}";
        var correlationLengthXString = $"XiX = {correlationLengthX} +- {correlationLengthXSigma}";
        var correlationLengthYString = $"XiY = {correlationLengthY} +- {correlationLengthYSigma}";
        var renormalisedCorrelationLengthString = $" Xi = {renormalisedCorrelationLength} +- {renormalisedCorrelationLengthSigma}";
        var susceptibilityString = $"Chi = {susceptibility} +- {susceptibilitySigma}";

        Console.WriteLine(magnetisationString);
        Console.WriteLine(magnetisationSquaredString);
        Console.WriteLine(magnetisationAbsoluteString + "\n");
        Console.WriteLine(energyString + "\n");
        Console.WriteLine(correlationLengthXString);
        Console.WriteLine(correlationLengthYString);
        Console.WriteLine(renormalisedCorrelationLengthString);
        Console.WriteLine(susceptibilityString);

        if (saveMeasurements)
        {
            writeMeasurements();
        }

        return new Observables(
            magnetisation,
            magnetisationSquared,
            magnetisationAbsolute,
            energy,
            susceptibility,
            correlationLengthX,
            correlationLengthY,
            renormalisedCorrelationLength,
            magnetisationSigma,
            magnetisationSquaredSigma,
            magnetisationAbsoluteSigma,
            energySigma,
            susceptibilitySigma,
            correlationLengthXSigma,
            correlationLengthYSigma,
            renormalisedCorrelationLengthSigma,
            magnetisationList,
            magnetisationSquaredList,
            magnetisationAbsoluteList,
            energyList,
            correlationLengthXList,
            correlationLengthYList,
            renormalisedCorrelationLengthList,
            susceptibilityList);

        void writeMeasurements()
        {
            var measurementDataDirectory = FileHelpers.GetDataRootDirectory(new string[]
            {
            Convert.ToString(_latticeLength),
            "measurements"
            });

            if (!Directory.Exists(measurementDataDirectory))
            {
                Directory.CreateDirectory(measurementDataDirectory);
            }

            var completePath =
                Path.GetFullPath(Path.Combine(measurementDataDirectory, $"{_temperature:0.00000}", ".txt"));

            var results = new List<string>
            {
                magnetisationString,
                magnetisationSquaredString,
                magnetisationAbsoluteString,
                energyString,
                correlationLengthXString,
                correlationLengthYString,
                renormalisedCorrelationLengthString,
                susceptibilityString,
                "\n",
                $"m = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", magnetisationList) + "}}",
                $"mSquared = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", magnetisationSquaredList) + "}}",
                $"mAbs = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", magnetisationAbsoluteList) + "}}",
                $"energy = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", energyList) + "}}",
                $"xi_x = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", correlationLengthXList) + "}}",
                $"xi_y = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", correlationLengthYList) + "}}",
                $"xi = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", renormalisedCorrelationLengthList) + "}}",
                $"chi = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", susceptibilityList) + "}}"
            };

            File.WriteAllLines(completePath, results);
        }
    }

    public void Thermalise(
        int thermalisationStepsInLatticeSizeUnit = ThermalisationStepsInMCSweepUnit,
        SpinUpdateMethod spinUpdateMethod = SpinUpdateMethod.Wolff,
        bool saveLattice = true,
        bool resetIterationCountDuringSave = false)
    {
        Simulation.RunMonteCarlo(
            _beta,
            _j,
            _h,
            thermalisationStepsInLatticeSizeUnit * Simulation.TotalSpinsCount,
            spinUpdateMethod,
            _jY,
            _randomSeed);

        if (saveLattice)
        {
            var iterationSteps = resetIterationCountDuringSave ? thermalisationStepsInLatticeSizeUnit
                                                                 * Simulation.TotalSpinsCount
                                                                : _previousIterationCount
                                                                  + thermalisationStepsInLatticeSizeUnit
                                                                  * Simulation.TotalSpinsCount;
            FileHelpers.SaveSpinConfiguration(
                Simulation.Lattice.Spins,
                _temperature,
                iterationSteps,
                isInByte: false);
        }

        Console.WriteLine(
            $"M (T={_temperature}) = {Simulation.Hamiltonian.GetAverageMagnetisation(_j, _h)}"
          + " after thermalisation.");
    }

    private static double GetVariance(List<double> measurements, double mean)
    {
        var measurementsCount = measurements.Count;

        return 1.0 / (measurementsCount - 1.0) * measurements.Select(m => (m - mean) * (m - mean)).Sum();
    }
}
