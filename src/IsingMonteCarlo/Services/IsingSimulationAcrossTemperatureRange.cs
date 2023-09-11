using System.Globalization;
using System.Transactions;

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
    }

    public void RunWithMeasurementsAcrossTemperatureRange(
        IEnumerable<double> temperatures,
        int iterationStepsBetweenMeasurements,
        int measurementsCount = 20,
        int measurementsRepetitionCount = 20,
        int thermalisationStepsInMCSweepUnit = 100_000,
        string? prethermalisedFirstLatticeFile = null,
        bool usePreviousTemperatureSpinsAsInitialConfiguration = true,
        bool loadFileForEachTemperature = false,
        bool saveLattice = true,
        bool saveMeasurements = true)
    {
        var enumeratedTemperatures = temperatures?.ToList() ?? throw new ArgumentNullException(nameof(temperatures));

        var iterationCount = 0;
        var previousTemperature = 0.0;
        var thermalisationMonteCarloSweepCount = thermalisationStepsInMCSweepUnit;

        var magnetisationList = new List<List<double>>();
        var magnetisationSquaredList = new List<List<double>>();
        var magnetisationAbsoluteList = new List<List<double>>();
        var energyList = new List<List<double>>();
        var correlationLengthXList = new List<List<double>>();
        var correlationLengthYList = new List<List<double>>();
        var renormalisedCorrelationLengthList = new List<List<double>>();
        var susceptibilityList = new List<List<double>>();

        foreach (var temperature in enumeratedTemperatures)
        {
            string? filename = null;
            if (loadFileForEachTemperature && prethermalisedFirstLatticeFile is null)
            {
                filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(
                    _latticeLength,
                    temperature);
            }
            else if (prethermalisedFirstLatticeFile is not null)
            {
                if (iterationCount is 0)
                {
                    filename = FileHelpers.GetDataRootDirectory(
                        new List<string> { _latticeLength.ToString(), prethermalisedFirstLatticeFile });
                    thermalisationStepsInMCSweepUnit = 0;
                }
                else
                {
                    filename = LoadExistingFileOrThermaliseFromPreviousTemperatureConfiguration(temperature, previousTemperature);
                }
            }

            var simulationWithObservables = new IsingSimulationWithObservablesComputation(
                filename,
                _dimension,
                _latticeLength,
                temperature,
                _j,
                _h,
                _spinUpdateMethod,
                _jY,
                randomSeed: _randomSeed);

            Console.Write($"\nStart run with T={temperature}...\n");
            Console.Write("\r");

            // iterationStepsBetweenMeasurements (suggested): _totalSpinsCount * [20, 100]
            var observables = simulationWithObservables.RunWithMeasurements(
                iterationStepsBetweenMeasurements,
                measurementsCount,
                measurementsRepetitionCount,
                thermalisationStepsInMCSweepUnit,
                saveLattice,
                saveMeasurements,
                resetIterationCountDuringSave: true);

            calculateExpectationValues(observables, temperature);

            Console.WriteLine($"\nRun with T={temperature} completed.\n");

            thermalisationStepsInMCSweepUnit = thermalisationMonteCarloSweepCount;
            previousTemperature = temperature;

            ++iterationCount;
        }

        writeMeasurements(enumeratedTemperatures);

        void calculateExpectationValues(Observables observables, double temperature)
        {
            magnetisationList.Add(
                new List<double>(3) { temperature, observables.Magnetisation, observables.MagnetisationSigma });
            magnetisationSquaredList.Add(
                new List<double>(3)
                {
                    temperature, observables.MagnetisationSquared, observables.MagnetisationSquaredSigma
                });
            magnetisationAbsoluteList.Add(
                new List<double>(3)
                {
                    temperature, observables.MagnetisationAbsolute, observables.MagnetisationAbsoluteSigma
                });
            energyList.Add(new List<double>(3) { temperature, observables.Energy, observables.EnergySigma });
            correlationLengthXList.Add(
                new List<double>(3)
                {
                    temperature,
                    observables.CorrelationLengthX,
                    observables.CorrelationLengthXSigma
                });
            correlationLengthYList.Add(
                new List<double>(3)
                {
                    temperature,
                    observables.CorrelationLengthY,
                    observables.CorrelationLengthYSigma
                });
            renormalisedCorrelationLengthList.Add(
                new List<double>(3)
                {
                    temperature,
                    observables.RenormalisedCorrelationLength,
                    observables.RenormalisedCorrelationLengthSigma
                });
            susceptibilityList.Add(
                new List<double>(3) { temperature, observables.Susceptibility, observables.SusceptibilitySigma });
        }

        void writeMeasurements(List<double> boltzmannTemperatures)
        {
            var measurementDataDirectory =
                FileHelpers.GetDataRootDirectory(new[] { Convert.ToString(_latticeLength), "measurements" });

            if (!Directory.Exists(measurementDataDirectory))
            {
                Directory.CreateDirectory(measurementDataDirectory);
            }

            var completePath = Path.GetFullPath(
                Path.Combine(
                    measurementDataDirectory,
                    $"T=[{boltzmannTemperatures.First():0.00000}, {boltzmannTemperatures.Last():0.00000}].txt"));

            var results = new List<string>
            {
                "m = "
              + "{" + string.Join(", ", magnetisationList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "mSquared = "
              + "{" + string.Join(", ", magnetisationSquaredList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "mAbs = "
              + "{" + string.Join(", ", magnetisationAbsoluteList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "energy = "
              + "{" + string.Join(", ", energyList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "xi_x = "
              + "{" + string.Join(", ", correlationLengthXList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "xi_y = "
              + "{" + string.Join(", ", correlationLengthYList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "xi = "
              + "{" + string.Join(", ", renormalisedCorrelationLengthList.Select(el => "{" + string.Join(", ", el) + "}")) + "}",
                "chi = "
              + "{" + string.Join(", ", susceptibilityList.Select(el => "{" + string.Join(", ", el) + "}")) + "}"
            };

            File.WriteAllLines(completePath, results);
        }

        string LoadExistingFileOrThermaliseFromPreviousTemperatureConfiguration(double temperature, double temperatureFromPreviousRun)
        {
            var currentTemperatureLatticeFilename = FileHelpers.GetFilename(
                _latticeLength,
                temperature,
                thermalisationStepsInMCSweepUnit);
            var dataDirectory = FileHelpers.GetDataLatticeLengthSubdirectory(_latticeLength);
            var currentTemperatureLatticeFile =
                Path.GetFullPath(Path.Combine(dataDirectory, currentTemperatureLatticeFilename, ".bin"));
            var previousTemperatureLatticeFilename = FileHelpers.GetFilename(
                _latticeLength,
                temperatureFromPreviousRun,
                thermalisationStepsInMCSweepUnit);
            var previousTemperatureLatticeFile =
                Path.GetFullPath(Path.Combine(dataDirectory, previousTemperatureLatticeFilename, ".bin"));

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
            var simulationWithObservables = new IsingSimulationWithObservablesComputation(
                filename: null,
                _dimension,
                _latticeLength,
                temperature,
                _j,
                _h,
                _spinUpdateMethod,
                _jY,
                randomSeed: _randomSeed);

            simulationWithObservables.Thermalise(
                thermalisationStepsInMCSweepUnit,
                _spinUpdateMethod,
                saveLattice: true);

            var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(simulationWithObservables.Simulation.Lattice.Spins);
            var resizedBitmap = DrawHelpers.ResizeBitmap(bitmap);
            DrawHelpers.SaveBitmapAsPNG(resizedBitmap, _latticeLength, temperature, resize: false);

            Console.WriteLine($"Run with T={temperature} completed.\n");
        }
    }
}
