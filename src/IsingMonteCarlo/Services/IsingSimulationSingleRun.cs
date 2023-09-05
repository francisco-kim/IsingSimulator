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

    public IsingSimulationSingleRun(
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

        Magnetisation = double.NaN;
        MagnetisationSquared = double.NaN;
        MagnetisationAbsolute = double.NaN;
        Energy = double.NaN;
        Susceptibility = double.NaN;
        RenormalisedCorrelationLength = double.NaN;

        MagnetisationSigma = double.NaN;
        MagnetisationSquaredSigma = double.NaN;
        MagnetisationAbsoluteSigma = double.NaN;
        EnergySigma = double.NaN;
        SusceptibilitySigma = double.NaN;
        RenormalisedCorrelationLengthSigma = double.NaN;

        MagnetisationList = new List<double>();
        MagnetisationSquaredList = new List<double>();
        MagnetisationAbsoluteList = new List<double>();
        EnergyList = new List<double>();
        SusceptibilityList = new List<double>();
        RenormalisedCorrelationLengthList = new List<double>();
        CorrelationLengthXList = new List<double>();
        CorrelationLengthYList = new List<double>();

        var (initialSpinConfiguration, _, previousIterationCount) =
            SpinConfigurationBuilder.InitialiseLattice(
                filename,
                dimension,
                latticeLength,
                initialSpinDownRatio,
                randomSeed);

        _previousIterationCount = previousIterationCount;

        Simulation = new IsingMonteCarloSimulation(dimension, latticeLength, initialSpinConfiguration);
    }

    public IsingMonteCarloSimulation Simulation { get; }

    public double Magnetisation { get; private set; }

    public double MagnetisationSquared { get; private set; }

    public double MagnetisationAbsolute { get; private set; }

    public double Energy { get; private set; }

    public double Susceptibility { get; private set; }

    public double RenormalisedCorrelationLength { get; private set; }

    public double MagnetisationSigma { get; private set; }

    public double MagnetisationSquaredSigma { get; private set; }

    public double MagnetisationAbsoluteSigma { get; private set; }

    public double EnergySigma { get; private set; }

    public double SusceptibilitySigma { get; private set; }

    public double RenormalisedCorrelationLengthSigma { get; private set; }

    public List<double> MagnetisationList { get; }

    public List<double> MagnetisationSquaredList { get; }

    public List<double> MagnetisationAbsoluteList { get; }

    public List<double> EnergyList { get; }

    public List<double> SusceptibilityList { get; }

    public List<double> RenormalisedCorrelationLengthList { get; }

    public List<double> CorrelationLengthXList { get; }

    public List<double> CorrelationLengthYList { get; }

    public void RunWithMeasurements(
        int iterationStepsBetweenMeasurements,
        int measurementsCount,
        int measurementsRepetitionCount,
        int thermalisationStepsInLatticeSizeUnit = ThermalisationStepsInLatticeSizeUnit,
        bool saveLattice = true,
        bool saveMeasurements = true,
        bool resetIterationCountDuringSave = false)
    {
        Thermalise(
            thermalisationStepsInLatticeSizeUnit,
            _spinUpdateMethod,
            saveLattice,
            resetIterationCountDuringSave);

        MeasurementsRun(
            iterationStepsBetweenMeasurements,
            measurementsCount,
            measurementsRepetitionCount,
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

    public void MeasurementsRun(
        int iterationStepsBetweenMeasurements,
        int measurementsCount,
        int measurementsRepetitionCount,
        bool saveMeasurements = true)
    {
        for (var measurementRepetition = 0;
             measurementRepetition < measurementsRepetitionCount;
             measurementRepetition++)
        {
            Simulation.RunMonteCarloWithObservablesComputation(
                _beta,
                _j,
                _h,
                iterationStepsBetweenMeasurements,
                measurementsCount,
                _spinUpdateMethod,
                _jY,
                _randomSeed);

            MagnetisationList.AddRange(Simulation.MagnetisationList);
            MagnetisationSquaredList.AddRange(Simulation.MagnetisationSquaredList);
            MagnetisationAbsoluteList.AddRange(Simulation.MagnetisationAbsoluteList);
            EnergyList.AddRange(Simulation.EnergyList);

            SusceptibilityList.Add(Simulation.Susceptibility);

            if (!double.IsNaN(Simulation.CorrelationLengthX))
            {
                CorrelationLengthXList.Add(Simulation.CorrelationLengthX);
            }

            if (!double.IsNaN(Simulation.CorrelationLengthY))
            {
                CorrelationLengthYList.Add(Simulation.CorrelationLengthY);
            }

            if (!double.IsNaN(Simulation.RenormalisedCorrelationLength))
            {
                RenormalisedCorrelationLengthList.Add(Simulation.RenormalisedCorrelationLength);
            }

            Simulation.ResetObservables();
        }

        Magnetisation = MagnetisationList.Sum() / MagnetisationList.Count;
        MagnetisationSigma = Math.Sqrt(GetVariance(MagnetisationList, Magnetisation));

        MagnetisationSquared = MagnetisationSquaredList.Sum() / MagnetisationSquaredList.Count;
        MagnetisationSquaredSigma = Math.Sqrt(GetVariance(MagnetisationSquaredList, MagnetisationSquared));

        MagnetisationAbsolute = MagnetisationAbsoluteList.Sum() / MagnetisationAbsoluteList.Count;
        MagnetisationAbsoluteSigma = Math.Sqrt(GetVariance(MagnetisationAbsoluteList, MagnetisationAbsolute));

        Energy = EnergyList.Sum() / EnergyList.Count;
        EnergySigma = Math.Sqrt(GetVariance(EnergyList, Energy));

        Console.WriteLine($" M  = {Magnetisation} +- {MagnetisationSigma}");
        Console.WriteLine($"M^2 = {MagnetisationSquared} +- {MagnetisationSquaredSigma}");
        Console.WriteLine($"|M| = {MagnetisationAbsolute} +- {MagnetisationAbsoluteSigma}");

        Susceptibility = SusceptibilityList.Sum() / SusceptibilityList.Count;
        SusceptibilitySigma = Math.Sqrt(GetVariance(SusceptibilityList, Susceptibility));

        RenormalisedCorrelationLength =
            RenormalisedCorrelationLengthList.Sum() / RenormalisedCorrelationLengthList.Count;
        RenormalisedCorrelationLengthSigma =
            Math.Sqrt(GetVariance(RenormalisedCorrelationLengthList, RenormalisedCorrelationLength));

        Console.WriteLine($"Chi = {Susceptibility} +- {SusceptibilitySigma}");
        Console.WriteLine($" Xi = {RenormalisedCorrelationLength} +- {RenormalisedCorrelationLengthSigma}");

        if (saveMeasurements)
        {
            SaveMeasurements();
        }
    }

    public void Thermalise(
        int thermalisationStepsInLatticeSizeUnit = ThermalisationStepsInLatticeSizeUnit,
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
            LatticeConfigurationSaver.SaveLattice(
                Simulation.Lattice.Spins,
                _temperature,
                iterationSteps,
                isBinary: false);
        }

        Console.WriteLine(
            $"M (T={_temperature}) = {Simulation.Hamiltonian.GetAverageMagnetisation(_j, _h)}"
          + " after thermalisation.\n");
    }

    private void SaveMeasurements()
    {
        var (_, measurementDataDirectory) = LatticeConfigurationSaver.GetFilename(
            Simulation.Lattice.Spins,
            _temperature,
            iterationCount: 0);
        measurementDataDirectory = Path.GetFullPath(Path.Combine(measurementDataDirectory, path2: "measurements"));

        if (!Directory.Exists(measurementDataDirectory))
        {
            Directory.CreateDirectory(measurementDataDirectory);
        }

        var completePathWithoutFileExtension =
            Path.GetFullPath(Path.Combine(measurementDataDirectory, $"{_temperature:0.0000}"));

        var results = new List<string>
        {
            $"m = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", MagnetisationList) + "}}",
            $"mSquared = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", MagnetisationSquaredList) + "}}",
            $"mAbs = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", MagnetisationAbsoluteList) + "}}",
            $"energy = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", EnergyList) + "}}",
            $"chi = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", SusceptibilityList) + "}}",
            $"xi = " + "{" + $"{_temperature}, " + "{" + string.Join(", ", RenormalisedCorrelationLengthList) + "}}"
        };

        File.WriteAllLines(
            completePathWithoutFileExtension + ".num",
            results);
    }

    private static double GetVariance(List<double> measurements, double mean)
    {
        var measurementsCount = measurements.Count;

        return 1.0 / (measurementsCount - 1.0) * measurements.Select(m => (m - mean) * (m - mean)).Sum();
    }
}
