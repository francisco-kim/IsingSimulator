using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Services;

//if (args.Length != 1)
//{
//    throw new ArgumentException(nameof(args));
//}

//var choice = Convert.ToInt32(args.First());
Console.WriteLine("0: single Monte Carlo run");
Console.WriteLine("1: Thermalisation (single run)");
Console.WriteLine("2: Monte Carlo run across temperature-range");
Console.WriteLine("3: Thermalisation across temperature-range");
Console.WriteLine("4: Generate png image");
Console.WriteLine("5: Renormalise png image");
Console.WriteLine("6: Continue thermalisation of 3^11\n");
Console.Write("Choice: ");

//var choice = 1;
var choice = Convert.ToInt32(Console.ReadLine());

var latticeLength = 19683;
var boltzmannTemperature = 2.2692;
var spinUpdateMethod = SpinUpdateMethod.Wolff;
if (choice is not 6)
{
    //var latticeLength = 19683;
    Console.Write("Lattice length: ");
    var latticeLengthInput = Console.ReadLine();
    // T = 2.2691853 = 2 / ln(1 + sqrt(2));
    latticeLength = latticeLengthInput is "" ? latticeLength : Convert.ToInt32(latticeLengthInput);
    Console.Write($"{latticeLength}\n");

    //var boltzmannTemperature = 2.0;
    Console.Write("Temperature: ");
    var temperatureInput = Console.ReadLine();
    // T = 2.2691853 = 2 / ln(1 + sqrt(2));
    boltzmannTemperature = temperatureInput is "" ? boltzmannTemperature : Convert.ToInt32(temperatureInput);
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

var iterationStepsBetweenMeasurements = totalSpinsCount * 5;
const int measurementsCount = 40;
const int measurementsRepetitionCount = 30;

if (choice == 0 || choice == 1)
{
    //string? filename = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";
    string? filename = null;

    var thermalisationStepsInLatticeSizeUnit =
        (filename is null) ? 100_000 : 0;

    var singleRunSimulation = new IsingSimulationSingleRun(filename,
                                                           dimension,
                                                           latticeLength,
                                                           boltzmannTemperature,
                                                           j,
                                                           h,
                                                           spinUpdateMethod,
                                                           randomSeed: randomSeed);

    if (choice == 0)
    {
        singleRunSimulation.RunWithMeasurements(
            iterationStepsBetweenMeasurements,
            measurementsCount,
            measurementsRepetitionCount,
            thermalisationStepsInLatticeSizeUnit,
            saveLattice: true,
            saveMeasurements: true);

        var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(singleRunSimulation.Simulation.Lattice.Spins);

        DrawHelpers.SaveBitmapAsPNG(bitmap, $"{latticeLength}_{boltzmannTemperature:0.0000}", resize: false);
    }
    else
    {
        singleRunSimulation.Thermalise(thermalisationStepsInLatticeSizeUnit, spinUpdateMethod, saveLattice: true);
    }
}
if (choice == 2 || choice == 3)
{
    var prethermalisedLatticeFile = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";
    var thermalisationStepsInLatticeSizeUnit = 100_000;

    var rangeCount = 8;
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

    if (choice == 2)
    {
        temperatureRangeSimulation.ThermaliseAcrossTemperatureRange(
            temperatures,
            thermalisationStepsInLatticeSizeUnit,
            spinUpdateMethod);

        Environment.Exit(0);
    }

    temperatureRangeSimulation.RunWithMeasurementsAcrossTemperatureRange(
        temperatures,
        iterationStepsBetweenMeasurements,
        measurementsCount,
        measurementsRepetitionCount,
        thermalisationStepsInLatticeSizeUnit,
        prethermalisedLatticeFile,
        usePreviousTemperatureSpinsAsInitialConfiguration: true,
        loadFile: true,
        saveLattice: true,
        saveMeasurements: true);
}
if (choice == 4)
{
    //var filename = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";
    var filename = $"81_2.2800_656100000.dat";

    var (initialSpinConfiguration, _, _) =
        SpinConfigurationBuilder.InitialiseLattice(
            filename,
            dimension,
            latticeLength);

    var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(initialSpinConfiguration);

    DrawHelpers.SaveBitmapAsPNG(bitmap, $"{latticeLength}_{boltzmannTemperature:0.0000}", resize: false);

    Environment.Exit(0);
}
if (choice == 5)
{
    var filename = $"81_{boltzmannTemperature:0.0000}_656100000.dat";
    var temperature = boltzmannTemperature;

    var initialLattice = FileHelpers.LoadLattice(filename, dimension);
    //Renormaliser.GenerateRenormalisedLatticeSpinConfigurationImages(initialLattice, twoPowNine, resize: true);
    Renormaliser.GenerateRenormalisedLatticeSpinConfigurationImages(initialLattice, 9, temperature, resize: true);
}
if (choice == 6)
{
    //var boltzmannTemperature = 2.0;
    Console.Write("100 MC-Sweep Unit: ");
    var thermalisationStepsInMillionMCSweepsUnitInput = Console.ReadLine();
    var thermalisationStepsIn100MCSweepUnit = thermalisationStepsInMillionMCSweepsUnitInput is ""
                                                       ? 1
                                                       : Convert.ToInt64(thermalisationStepsInMillionMCSweepsUnitInput);
    Console.Write($"{thermalisationStepsIn100MCSweepUnit}");

    string? filename = null;
    try
    {
        filename = FileHelpers.GetFirstDataFileWithLatticeSizeAndTemperature(latticeLength, boltzmannTemperature);
    }
    catch (Exception e)
    {
        // ignored
    }

    //string? filename = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";

    var singleRunSimulation = new IsingSimulationSingleRun(filename,
                                                           dimension,
                                                           latticeLength,
                                                           boltzmannTemperature,
                                                           j,
                                                           h,
                                                           spinUpdateMethod,
                                                           randomSeed: randomSeed);

    singleRunSimulation.ThermaliseLargeLattice(thermalisationStepsIn100MCSweepUnit, spinUpdateMethod: spinUpdateMethod);

    var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(singleRunSimulation.Simulation.Lattice.Spins);
    DrawHelpers.SaveBitmapAsPNG(bitmap, $"{latticeLength}_{boltzmannTemperature:0.0000}", resize: false);
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
