using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;
using IsingMonteCarlo.Services;

//if (args.Length != 1)
//{
//    throw new ArgumentException(nameof(args));
//}

//var choice = Convert.ToInt32(args.First());
Console.WriteLine("0: single Monte Carlo run");
Console.WriteLine("1: Thermalisation (single run)");
Console.WriteLine("2: Thermalisation across temperature-range (2D with J = -1)");
Console.WriteLine("3: Monte Carlo run across temperature-range");
Console.WriteLine("4: Thermalisation across temperature-range");
Console.WriteLine("5: Generate png image");
Console.WriteLine("6: Renormalise png image");
Console.WriteLine("7: Continue thermalisation of 3^11\n");
Console.Write("Choice: ");

//var choice = 1;
var choice = Convert.ToInt32(Console.ReadLine());

var latticeLength = 19683;
var boltzmannTemperature = 2.26923;    // T_c = 2.26919    T_c(L = 128) = 2.27557
var spinUpdateMethod = SpinUpdateMethod.Wolff;
if (choice is not 7)
{
    //var latticeLength = 19683;
    Console.Write("Lattice length: ");
    var latticeLengthInput = Console.ReadLine();
    // T = 2.2691853 = 2 / ln(1 + sqrt(2));
    latticeLength = latticeLengthInput is "" ? latticeLength : Convert.ToInt32(latticeLengthInput);
    Console.Write($"{latticeLength}\n");
    //var boltzmannTemperature = 2.0;
    Console.Write("Temperature (T_c = 2.26919; T_c(L = 128) = 2.27557): ");
    var temperatureInput = Console.ReadLine();
    // T = 2.2691853 = 2 / ln(1 + sqrt(2));
    boltzmannTemperature = temperatureInput is "" ? boltzmannTemperature : Convert.ToDouble(temperatureInput);
    Console.Write($"{boltzmannTemperature}\n");

    //var spinUpdateMethod = SpinUpdateMethod.Wolff;
    Console.Write("(g)lauber or (w)olff: ");
    var spinUpdateMethodInput = Console.ReadLine();
    spinUpdateMethod = spinUpdateMethodInput is "w" ? spinUpdateMethod : SpinUpdateMethod.Glauber;
    Console.Write($"{spinUpdateMethodInput}\n");
}

var dimension = 2;
var totalSpinsCount = latticeLength * latticeLength;
var j = -1.0;
var h = 0.0;
int? randomSeed = 41;
//var spinUpdateMethod = SpinUpdateMethod.Glauber;
//var spinUpdateMethod = SpinUpdateMethod.Wolff;

var iterationNeededForSingleChiXiMeasurement = totalSpinsCount / 2;
const int measurementsCountForChiXiExpectationValue = 200;
const int measurementsRepetitionCountForChiXiVariance = 5;

if (choice == 0 || choice == 1)
{
    string? filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
    // string? filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_1000000000.dat";
    // string? filename = null;

    var thermalisationStepsInLatticeSizeUnit =
        (filename is null) ? 100_000 : 0;

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
}
if (choice == 2)
{
    var temperatures = new List<double>()
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
    //var rangeCount = 8;
    //var deltaTemperature = 0.005;
    //var temperatures = Enumerable.Range(0, rangeCount)
    //                             .Select(i => Math.Round((0.005 * i + 2.26) / deltaTemperature) * deltaTemperature);

    //string? filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_1000000000.dat";
    //try
    //{
    //    filename = FileHelpers.GetFirstDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
    //}
    //catch (Exception)
    //{
    //    // ignored
    //}
    foreach (var temperature in temperatures)
    {
        string? filename = null;

        var thermalisedSpins = Fast2DIsingMonteCarloSimulator.ThermaliseLargeLattice(
            filename,
            latticeLength,
            temperature,
            thermalisationStepsIn100MCSweepUnit: 1000,
            spinUpdateMethod,
            resetIterationCountDuringSave: false,
            randomSeed: null,
            verbose: false);

        var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(thermalisedSpins);
        DrawHelpers.SaveBitmapAsPNG(bitmap, latticeLength, boltzmannTemperature, resize: true);

        Console.WriteLine($"T = {temperature} completed.");
    }
}
if (choice == 3 || choice == 4)
{
    // var prethermalisedLatticeFile = $"{latticeLength}_{boltzmannTemperature:0.00000}_1000.dat";
    var prethermalisedLatticeFile = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature)
                                    ?? "";
    var thermalisationStepsInLatticeSizeUnit = 100_000;

    //var rangeCount = 8;
    //var deltaTemperature = 0.005;
    //var temperatures = Enumerable.Range(0, rangeCount)
    //                             .Select(i => Math.Round((0.005 * i + 2.26) / deltaTemperature) * deltaTemperature);
    var temperatures = new List<double>()
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
            temperatures,
            thermalisationStepsInLatticeSizeUnit);

        Environment.Exit(0);
    }

    temperatureRangeSimulation.RunWithMeasurementsAcrossTemperatureRange(
        temperatures,
        iterationNeededForSingleChiXiMeasurement,
        measurementsCountForChiXiExpectationValue,
        measurementsRepetitionCountForChiXiVariance,
        thermalisationStepsInLatticeSizeUnit,
        prethermalisedLatticeFile,
        usePreviousTemperatureSpinsAsInitialConfiguration: false,
        loadFileForEachTemperature: true,
        saveLattice: true,
        saveMeasurements: true);
}
if (choice == 5)
{
    // var filename = $"81_2.2800_656100000.dat";
    // var filename = $"19683_2.26920_1.bin";
    var filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);

    var (initialSpinConfiguration, _, _) =
        SpinConfigurationBuilder.InitialiseLattice(
            filename,
            dimension,
            latticeLength);

    var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(initialSpinConfiguration);

    DrawHelpers.SaveBitmapAsPNG(bitmap, latticeLength, boltzmannTemperature, resize: true);

    Environment.Exit(0);
}
if (choice == 6)
{
    var filesInLatticeLengthFolder = Directory.GetFiles(FileHelpers.GetDataLatticeLengthSubdirectory(latticeLength))
        .Where(file => file.Split('.').Last() == "dat"
            || file.Split('.').Last() == "bin");
    foreach (string file in filesInLatticeLengthFolder)
    {
        Console.WriteLine(file);
    }
    Console.WriteLine("Enter file name (.dat or .bin): ");

    var filename = Console.ReadLine() ?? throw new ArgumentNullException();
    // var filename = $"243_2.26920_1609932704.dat";
    // var filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_656100000.dat";

    var initialLattice = FileHelpers.LoadLattice(filename, dimension);
    //Renormaliser.GenerateRenormalisedLatticeSpinConfigurationImages(initialLattice, twoPowNine, resize: true);
    Renormaliser.GenerateRenormalisedLatticeSpinConfigurationImages(initialLattice, 9, boltzmannTemperature, resize: true);
}
if (choice == 7)
{
    Console.Write("100 MC-Sweep Unit: ");
    var thermalisationStepsInMillionMCSweepsUnitInput = Console.ReadLine();
    var thermalisationStepsIn100MCSweepUnit = thermalisationStepsInMillionMCSweepsUnitInput is ""
                                                       ? 1
                                                       : Convert.ToInt64(thermalisationStepsInMillionMCSweepsUnitInput);
    Console.Write($"{thermalisationStepsIn100MCSweepUnit}\n");

    // var filename = $"{latticeLength}_{boltzmannTemperature:0.00000}_1.bin";
    //string? filename = null;
    var filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
    try
    {
        filename = FileHelpers.GetLastDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
    }
    catch (Exception)
    {
        // ignored
    }

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
}

// else
// {
//     var criticalRangeCount = 8;

//     var temperaturesList = new List<double>(3) { 0.0001, 0.8, 1.6 };
//     var criticalTemperatureRange = Enumerable.Range(0, criticalRangeCount)
//                                              .Select(i => 0.5 / criticalRangeCount * i + 2.0);
//     var highTemperatureRange = new List<double>(5) { 2.5, 3.0, 3.5, 4.0, 4.5 };

//     temperaturesList.AddRange(criticalTemperatureRange);
//     temperaturesList.AddRange(highTemperatureRange);

//     var magnetisations = temperaturesList.Select((t, index)
//                                                  => getMagnetisationFromEachRun(1.0 / t, index));

//     var magnetisationData = temperaturesList.Zip(magnetisations);
//     var magnetisationDataList = magnetisationData.Select(data => $"{{{data.First}, {data.Second}}}");
//     Console.WriteLine(string.Join(", ", magnetisationDataList));

//     double getMagnetisationFromEachRun(double beta, int iterationIndex)
//     {
//         monteCarlo.RunMonteCarlo(beta,
//                                  j,
//                                  h,
//                                  thermalisationCount,
//                                  spinUpdateMethod,
//                                  randomSeed: 17);

//         Console.WriteLine($"Iteration {iterationIndex} completed.");

//         return monteCarlo.Hamiltonian.GetAverageMagnetisation(j, h);
//     }
// }
