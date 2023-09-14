using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;
using IsingMonteCarlo.Services;

//if (args.Length != 1)
//{
//    throw new ArgumentException(nameof(args));
//}

//var choice = Convert.ToInt32(args.First());

//var choice = 1;
Console.WriteLine(value: "Ising Monte Carlo Simulator\n");

Console.WriteLine(value: "Option 0: Single Monte Carlo run");
Console.WriteLine(value: "Option 1: Thermalisation (single run)");
Console.WriteLine(value: "Option 2: Thermalisation across temperature-range (2D with J = -1)");
Console.WriteLine(value: "Option 3: Monte Carlo run across temperature-range");
Console.WriteLine(value: "Option 4: Thermalisation across temperature-range");
Console.WriteLine(value: "Option 5: Generate png image");
Console.WriteLine(value: "Option 6: Renormalise png image");
Console.WriteLine(value: "Option 7: Zoom-in with png image");
Console.WriteLine(value: "Option 8: Continue thermalisation of 3^11\n");

Console.Write(value: "Please choose one of the above options: ");
var choice = Convert.ToInt32(Console.ReadLine());
Console.WriteLine("\n");

var latticeLength = 128;
var boltzmannTemperature = 2.27557; // T_c = 2.2691853 = 2 / ln(1 + sqrt(2))   T_c(L = 128) = 2.27557
var spinUpdateMethod = SpinUpdateMethod.Wolff;

if (choice is not 6 && choice is not 7 && choice is not 8)
{
    Console.Write(value: "* Lattice length (max. 19683): ");
    var latticeLengthInput = Console.ReadLine();
    latticeLength = latticeLengthInput is "" ? latticeLength : Convert.ToInt32(latticeLengthInput);
    Console.Write($"      {latticeLength}\n");
}

if (choice is 0 or 1 or 5 or 8)
{
    Console.Write(value: "* Temperature (T_c = 2.26919 / T_c(L = 128) = 2.27557 / T_c(L = 19683) = 2.26923): ");
    var temperatureInput = Console.ReadLine();
    boltzmannTemperature = temperatureInput is "" ? boltzmannTemperature : Convert.ToDouble(temperatureInput);
    Console.Write($"      {boltzmannTemperature}\n");
}

if (choice is 0 or 1 or 2 or 3 or 4)
{
    Console.Write(value: "* (g)lauber or (w)olff: ");
    var spinUpdateMethodInput = Console.ReadLine();
    spinUpdateMethod = spinUpdateMethodInput is "g" ? SpinUpdateMethod.Glauber : SpinUpdateMethod.Wolff;
    Console.Write($"      Spin dynamics: {spinUpdateMethod}\n");
}

var dimension = 2;
var totalSpinsCount = latticeLength * latticeLength;
var j = -1.0;
var h = 0.0;
int? randomSeed = 41;

var iterationNeededForSingleChiXiMeasurement = 10 * totalSpinsCount;
const int measurementsCountForChiXiExpectationValue = 200;
const int measurementsRepetitionCountForChiXiVariance = 10;

switch (choice)
{
    case 0 or 1:
    {
        // string? filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_1000000000.dat";
        // string? filename = null;
        var filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
        if (filename is not null)
        {
            Console.WriteLine("Loaded " + filename + ".\n");
        }

        var thermalisationStepsInLatticeSizeUnit =
            filename is null ? 100_000 : 0;

        var simulationWithObservables = new IsingSimulationWithObservablesComputation(
            filename,
            dimension,
            latticeLength,
            boltzmannTemperature,
            j,
            h,
            spinUpdateMethod,
            randomSeed: randomSeed);

        if (choice == 0)
        {
            var _ = simulationWithObservables.RunWithMeasurements(
                iterationNeededForSingleChiXiMeasurement,
                measurementsCountForChiXiExpectationValue,
                measurementsRepetitionCountForChiXiVariance,
                thermalisationStepsInLatticeSizeUnit,
                saveLattice: true,
                saveMeasurements: true);

            var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(simulationWithObservables.Simulation.Lattice.Spins);

            DrawHelpers.SaveBitmapAsPNG(bitmap, latticeLength, boltzmannTemperature, resize: false);
        }
        else
        {
            simulationWithObservables.Thermalise(thermalisationStepsInLatticeSizeUnit, spinUpdateMethod, saveLattice: true);
        }

        break;
    }
    case 2:
    {
        var temperatures = new List<double>
        {
            2.26,
            2.265,
            2.27,
            2.275,
            2.28,
            2.285,
            2.29,
            2.295
        };

        foreach (var temperature in temperatures)
        {
            string? filename = null;

            var thermalisedSpins = Fast2DIsingMonteCarloSimulator.ThermaliseLargeLattice(
                filename,
                latticeLength,
                temperature,
                thermalisationStepsIn100MCSweepUnit: 1,
                spinUpdateMethod,
                resetIterationCountDuringSave: false,
                randomSeed: null,
                verbose: false);

            var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(thermalisedSpins);
            DrawHelpers.SaveBitmapAsPNG(bitmap, latticeLength, boltzmannTemperature, resize: true);

            Console.WriteLine($"T = {temperature} completed.");
        }

        break;
    }
    case 3 or 4:
    {
        Console.Write(value: "* Temperatures in the critical region? ((y)es or other): ");
        var isCriticalRegionInput = Console.ReadLine();
        var isCriticalRegion = isCriticalRegionInput is "y";
        Console.Write($"      {isCriticalRegion}\n");

        var thermalisationStepsInMCSweepUnit = 100_000;
        Console.Write(value: "* Thermalisation steps in Monte-Carlo sweep unit: ");
        var thermalisationStepsInMCSweepUnitInput = Console.ReadLine();
        thermalisationStepsInMCSweepUnit = thermalisationStepsInMCSweepUnitInput is ""
                                               ? thermalisationStepsInMCSweepUnit
                                               : Convert.ToInt32(thermalisationStepsInMCSweepUnitInput);
        if (thermalisationStepsInMCSweepUnitInput is "")
        {
            Console.Write($"      {thermalisationStepsInMCSweepUnit}\n");
        }

            // var prethermalisedLatticeFile = $"{latticeLength}_{boltzmannTemperature:0.00000}_1000.dat";
            //var prethermalisedLatticeFile = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
            string? prethermalisedLatticeFile = null;

        //var rangeCount = 8;
        //var deltaTemperature = 0.005;
        //var temperatures = Enumerable.Range(0, rangeCount)
        //                             .Select(i => Math.Round((0.005 * i + 2.26) / deltaTemperature) * deltaTemperature);
        var temperaturesInCriticalRegion = new List<double>
        {
            2.26,
            2.265,
            2.27,
            2.275,
            2.28,
            2.285,
            2.29,
            2.295
        };

        var temperatures = new List<double>
        {
            0.1,
            0.5,
            1.0,
            2.0,
            2.26,
            2.27,
            2.28,
            2.29,
            3.0,
            3.5,
            4.0
        };

        var boltzmannTemperatures = isCriticalRegion ? temperaturesInCriticalRegion : temperatures;

        var temperatureRangeSimulation = new IsingSimulationAcrossTemperatureRange(
            dimension,
            latticeLength,
            j,
            h,
            spinUpdateMethod,
            randomSeed: randomSeed);

        if (choice == 4)
        {
            temperatureRangeSimulation.ThermaliseAcrossTemperatureRange(
                boltzmannTemperatures,
                thermalisationStepsInMCSweepUnit);

            Environment.Exit(exitCode: 0);
        }

        temperatureRangeSimulation.RunWithMeasurementsAcrossTemperatureRange(
            boltzmannTemperatures,
            iterationNeededForSingleChiXiMeasurement,
            measurementsCountForChiXiExpectationValue,
            measurementsRepetitionCountForChiXiVariance,
            thermalisationStepsInMCSweepUnit,
            prethermalisedLatticeFile,
            usePreviousTemperatureSpinsAsInitialConfiguration: false,
            loadFileForEachTemperature: true,
            saveLattice: true,
            saveMeasurements: true);
        break;
    }
    case 5:
    {
        // var filename = $"81_2.2800_656100000.dat";
        // var filename = $"19683_2.26920_1.bin";
        var filename = getFilename(latticeLength, boltzmannTemperature);

        var (initialSpinConfiguration, _, _) =
            SpinConfigurationBuilder.InitialiseLattice(
                filename,
                dimension,
                latticeLength);

        var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(initialSpinConfiguration);

        DrawHelpers.SaveBitmapAsPNG(bitmap, latticeLength, boltzmannTemperature, resize: true);

        Environment.Exit(exitCode: 0);
        break;
    }
    case 6:
    {
        latticeLength = 19683;
        boltzmannTemperature = 2.26923;

        var filename = getFilename(latticeLength, boltzmannTemperature);

        // var filename = $"243_2.26920_1609932704.dat";
        // var filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_656100000.dat";

        var initialSpinConfiguration = FileHelpers.LoadSpinConfiguration(filename, out _, out _);
        Renormaliser.GenerateRenormalisedLatticeSpinConfigurationImages(
            initialSpinConfiguration,
            finalLatticeSizeLimit: 9,
            boltzmannTemperature,
            resize: true);
        break;
    }
    case 7:
    {
        latticeLength = 19683;
        boltzmannTemperature = 2.26923;

        var filename = getFilename(latticeLength, boltzmannTemperature);

        // var filename = $"243_2.26920_1609932704.dat";
        // var filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_656100000.dat";

        var initialSpinConfiguration = FileHelpers.LoadSpinConfiguration(filename, out _, out _);
        Renormaliser.GenerateZoomedLatticeSpinConfigurationImages(
            initialSpinConfiguration,
            finalLatticeSizeLimit: 9,
            boltzmannTemperature,
            resize: true);
        break;
    }
    case 8:
    {
        latticeLength = 19683;
        Console.Write($"* Lattice length: {latticeLength}\n");
        Console.Write($"* Temperature: {boltzmannTemperature}\n");
        Console.Write("* 100 MC-Sweep Unit: ");
        var thermalisationStepsIn100MCSweepsUnitInput = Console.ReadLine();
        Console.Write("\r");

        var thermalisationStepsIn100MCSweepUnit = 1.0;
        if (thermalisationStepsIn100MCSweepsUnitInput is "")
        {
            Console.Write($"      {thermalisationStepsIn100MCSweepUnit}\n");
        }
        else
        {
            thermalisationStepsIn100MCSweepUnit = Convert.ToDouble(thermalisationStepsIn100MCSweepsUnitInput);
        }

        var filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
        if (filename is not null)
        {
            Console.WriteLine($"Loaded {filename}.\n");
        }
        // filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_1.bin";

        var thermalisedSpins = Fast2DIsingMonteCarloSimulator.ThermaliseLargeLattice(
            filename,
            latticeLength,
            boltzmannTemperature,
            thermalisationStepsIn100MCSweepUnit,
            spinUpdateMethod,
            resetIterationCountDuringSave: false,
            randomSeed: null,
            verbose: true);

        var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(thermalisedSpins);
        DrawHelpers.SaveBitmapAsPNG(bitmap, latticeLength, boltzmannTemperature, resize: true);
        break;
    }
}

static string getFilename(int givenLatticeLength, double boltzmannTemperature)
{
    var filenameCandidate = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(givenLatticeLength, boltzmannTemperature);
    if (filenameCandidate is not null)
    {
        Console.Write($"Found {filenameCandidate}.\n" + "Proceed (enter) or load another file (any other character)?\n");
        var proceedConfirmation = Console.ReadKey().Key;

        if (proceedConfirmation is not ConsoleKey.Enter)
        {
            filenameCandidate = getFilenameFromUserInput(givenLatticeLength);
        }
    }
    else
    {
        filenameCandidate = getFilenameFromUserInput(givenLatticeLength);
    }

    return filenameCandidate;

    static string getFilenameFromUserInput(int latticeLengthToMatchFolderName)
    {
        var foundFile = false;
        var foundFilename = "";
        while (!foundFile)
        {
            Console.WriteLine(value: "\nExisting files: ");
            foundFile = FileHelpers.FoundFileInLatticeLengthFolder(
                latticeLengthToMatchFolderName,
                out var foundFilenameCandidate);
            foundFilename = foundFilenameCandidate ?? throw new ArgumentNullException(nameof(foundFilenameCandidate));
        }

        return foundFilename;
    }
}
