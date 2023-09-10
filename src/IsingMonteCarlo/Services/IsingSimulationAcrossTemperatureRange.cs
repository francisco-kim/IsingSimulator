using System.Globalization;

using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;

namespace IsingMonteCarlo.Services;

public sealed class IsingSimulationAcrossTemperatureRange
{
    private const int DefaultThermalisationStepsInMCSweepUnit = 20_000;

    private readonly SpinUpdateMethod _spinUpdateMethod;
    private readonly double _j;
    private readonly double _h;
    private readonly double? _jY;
    private readonly int _dimension;
    private readonly int _latticeLength;
    private readonly int? _randomSeed;

    public IsingSimulationAcrossTemperatureRange(
        int dimension,
        int latticeLength,
        double j,
        double h,
        SpinUpdateMethod spinUpdateMethod = SpinUpdateMethod.Wolff,
        double? jY = null,
        double initialSpinDownRatio = 0.25,
        int? randomSeed = 41)
    {
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

        MagnetisationList = new List<List<double>>();
        MagnetisationSquaredList = new List<List<double>>();
        MagnetisationAbsoluteList = new List<List<double>>();
        EnergyList = new List<List<double>>();
        SusceptibilityList = new List<List<double>>();
        RenormalisedCorrelationLengthList = new List<List<double>>();
    }

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

    public List<List<double>> MagnetisationList { get; private set; }

    public List<List<double>> MagnetisationSquaredList { get; private set; }

    public List<List<double>> MagnetisationAbsoluteList { get; private set; }

    public List<List<double>> EnergyList { get; private set; }

    public List<List<double>> SusceptibilityList { get; private set; }

    public List<List<double>> RenormalisedCorrelationLengthList { get; private set; }

    public void RunWithMeasurementsAcrossTemperatureRange(
        IEnumerable<double> temperatures,
        int iterationStepsBetweenMeasurements,
        int measurementsCount = 20,
        int measurementsRepetitionCount = 20,
        int thermalisationStepsInMCSweepUnit = 100_000,
        string prethermalisedFirstLatticeFile = "",
        bool usePreviousTemperatureSpinsAsInitialConfiguration = true,
        bool loadFileForEachTemperature = false,
        bool saveLattice = true,
        bool saveMeasurements = true)
    {
        var enumeratedTemperatures = temperatures?.ToList() ?? throw new ArgumentNullException(nameof(temperatures));

        var iterationCount = 0;
        var previousTemperature = 0.0;
        var thermalisationMonteCarloSweepCount = thermalisationStepsInMCSweepUnit;

        foreach (var temperature in enumeratedTemperatures)
        {
            string? filename = null;
            if (loadFileForEachTemperature && prethermalisedFirstLatticeFile is "")
            {
                filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(
                    _latticeLength,
                    temperature);
            }
            else if (prethermalisedFirstLatticeFile is not "")
            {
                var dataDirectory = FileHelpers.GetDataRootDirectory();

                if (iterationCount is 0)
                {
                    filename = Path.GetFullPath(Path.Combine(dataDirectory, prethermalisedFirstLatticeFile));
                    thermalisationStepsInMCSweepUnit = 0;
                }
                else
                {
                    filename = LoadExistingFileOrThermaliseFromPreviousTemperatureConfiguration(
                        dataDirectory,
                        temperature);
                }
            }

            var singleRunSimulation = new IsingSimulationSingleRun(
                filename,
                _dimension,
                _latticeLength,
                temperature,
                _j,
                _h,
                _spinUpdateMethod,
                _jY,
                randomSeed: _randomSeed);

            // iterationStepsBetweenMeasurements (suggested): _totalSpinsCount * [20, 100]
            singleRunSimulation.RunWithMeasurements(
                iterationStepsBetweenMeasurements,
                measurementsCount,
                measurementsRepetitionCount,
                thermalisationStepsInMCSweepUnit,
                saveLattice,
                saveMeasurements,
                resetIterationCountDuringSave: true);

            CalculateExpectationValues(singleRunSimulation, temperature);

            ResetObservables();

            Console.WriteLine($"\nRun with T={temperature} completed.\n");

            thermalisationStepsInMCSweepUnit = thermalisationMonteCarloSweepCount;
            previousTemperature = temperature;

            ++iterationCount;
        }

        SaveMeasurements(enumeratedTemperatures);

        string LoadExistingFileOrThermaliseFromPreviousTemperatureConfiguration(
            string dataDirectory,
            double temperature)
        {
            var thermalisationIterationCount =
                thermalisationStepsInMCSweepUnit * _latticeLength * _latticeLength;
            var currentTemperatureLatticeFile =
                Path.GetFullPath(
                    Path.Combine(
                        dataDirectory,
                        $"{_latticeLength}_{temperature:0.00000}_{thermalisationIterationCount.ToString(format: "G10", CultureInfo.InvariantCulture)}.dat"));
            var previousTemperatureLatticeFile =
                Path.GetFullPath(
                    Path.Combine(
                        dataDirectory,
                        $"{_latticeLength}_{previousTemperature:0.00000}_{thermalisationIterationCount.ToString(format: "G10", CultureInfo.InvariantCulture)}.dat"));

            var filenameToLoad = "";
            if (File.Exists(currentTemperatureLatticeFile))
            {
                filenameToLoad = currentTemperatureLatticeFile;
                thermalisationStepsInMCSweepUnit = 0;
            }

            if (File.Exists(previousTemperatureLatticeFile)
             && usePreviousTemperatureSpinsAsInitialConfiguration)
            {
                filenameToLoad = previousTemperatureLatticeFile;
            }

            return filenameToLoad;
        }
    }

    public void ThermaliseAcrossTemperatureRange(
        IEnumerable<double> temperatures,
        int thermalisationStepsInMCSweepUnit = DefaultThermalisationStepsInMCSweepUnit)
    {
        var enumeratedTemperatures = temperatures?.ToList() ?? throw new ArgumentNullException(nameof(temperatures));

        foreach (var temperature in enumeratedTemperatures)
        {
            var singleRunSimulation = new IsingSimulationSingleRun(
                filename: null,
                _dimension,
                _latticeLength,
                temperature,
                _j,
                _h,
                _spinUpdateMethod,
                _jY,
                randomSeed: _randomSeed);

            singleRunSimulation.Thermalise(
                thermalisationStepsInMCSweepUnit,
                _spinUpdateMethod,
                saveLattice: true);

            var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(singleRunSimulation.Simulation.Lattice.Spins);
            var resizedBitmap = DrawHelpers.ResizeBitmap(bitmap);
            DrawHelpers.SaveBitmapAsPNG(resizedBitmap, _latticeLength, temperature, resize: false);

            Console.WriteLine($"Run with T={temperature} completed.\n");
        }
    }

    private void CalculateExpectationValues(IsingSimulationSingleRun singleRunSimulation, double temperature)
    {
        MagnetisationList.Add(
            new List<double>(3)
            {
                temperature, singleRunSimulation.Magnetisation, singleRunSimulation.MagnetisationSigma
            });
        MagnetisationSquaredList.Add(
            new List<double>(3)
            {
                temperature, singleRunSimulation.MagnetisationSquared, singleRunSimulation.MagnetisationSquaredSigma
            });
        MagnetisationAbsoluteList.Add(
            new List<double>(3)
            {
                temperature,
                singleRunSimulation.MagnetisationAbsolute,
                singleRunSimulation.MagnetisationAbsoluteSigma
            });
        EnergyList.Add(
            new List<double>(3) { temperature, singleRunSimulation.Energy, singleRunSimulation.EnergySigma });
        SusceptibilityList.Add(
            new List<double>(3)
            {
                temperature, singleRunSimulation.Susceptibility, singleRunSimulation.SusceptibilitySigma
            });
        RenormalisedCorrelationLengthList.Add(
            new List<double>(3)
            {
                temperature,
                singleRunSimulation.RenormalisedCorrelationLength,
                singleRunSimulation.RenormalisedCorrelationLengthSigma
            });
    }

    private void SaveMeasurements(List<double> temperatures)
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
            Path.GetFullPath(
                Path.Combine(
                    measurementDataDirectory,
                    $"T=[{temperatures.First():0.00000}, {temperatures.Last():0.00000}]", ".txt"));

        var results = new List<string>
        {
            "m = " + "{" + string.Join(", ", MagnetisationList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
            "mSquared = " + "{" + string.Join(", ", MagnetisationSquaredList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
            "mAbs = " + "{" + string.Join(", ", MagnetisationAbsoluteList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
            "energy = " + "{" + string.Join(", ", EnergyList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
            "chi = " + "{" + string.Join(", ", SusceptibilityList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
            "xi = " + "{" + string.Join(", ", RenormalisedCorrelationLengthList.Select(el => "{" + string.Join(", ", el) + "}")) + "}"
        };

        File.WriteAllLines(completePath + ".num", results);
    }

    private void ResetObservables()
    {
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
    }

    private static double GetVariance(List<double> measurements, double mean)
    {
        var measurementsCount = measurements.Count;

        return 1.0 / (measurementsCount - 1.0) * measurements.Select(m => (m - mean) * (m - mean)).Sum();
    }
}
