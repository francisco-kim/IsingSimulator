using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Services;

var generatePNG = true;
var isSingleRun = false;
var thermaliseAcrossTemperatureRange = true;

var dimension = 2;
var latticeLength = 100;
var totalSpinsCount = latticeLength * latticeLength;
var j = -1.0;
var h = 0.0;
var spinUpdateMethod = SpinUpdateMethod.Glauber;
//var spinUpdateMethod = SpinUpdateMethod.Wolff;
int? randomSeed = 41;

var boltzmannTemperature = 2.28;

var iterationStepsBetweenMeasurements = totalSpinsCount * 20;
const int measurementsCount = 20;
const int measurementsRepetitionCount = 20;

if (generatePNG)
{
    string? filename = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";

    var (initialSpinConfiguration, _, _) =
        SpinConfigurationBuilder.InitialiseLattice(
            filename,
            dimension,
            latticeLength);

    var bitmap = DrawHelper.FromTwoDimIntArrayGray(initialSpinConfiguration);
    var resizedBitmap = DrawHelper.ResizeToLargerBitmap(bitmap, 512, 512);
    DrawHelper.SaveBmpAsPNG(resizedBitmap, $"{latticeLength}_{boltzmannTemperature:0.0000}");

    Environment.Exit(0);
}

if (isSingleRun)
{
    //string? filename = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";
    string? filename = null;

    var thermalisationStepsInLatticeSizeUnit = 0;
    if (filename is null)
    {
        thermalisationStepsInLatticeSizeUnit = 100_000;
    }

    var singleRunSimulation = new IsingSimulationSingleRun(filename,
                                                           dimension,
                                                           latticeLength,
                                                           boltzmannTemperature,
                                                           j,
                                                           h,
                                                           spinUpdateMethod,
                                                           randomSeed: randomSeed);

    singleRunSimulation.RunWithMeasurements(iterationStepsBetweenMeasurements,
                                     measurementsCount,
                                     measurementsRepetitionCount,
                                     thermalisationStepsInLatticeSizeUnit,
                                     saveLattice: true,
                                     saveMeasurements: true);

    var bitmap = DrawHelper.FromTwoDimIntArrayGray(singleRunSimulation.Simulation.Lattice.Spins);
    var resizedBitmap = DrawHelper.ResizeToLargerBitmap(bitmap, 512, 512);
    DrawHelper.SaveBmpAsPNG(resizedBitmap, $"{latticeLength}_{boltzmannTemperature:0.0000}");
}
else
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

    if (thermaliseAcrossTemperatureRange)
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