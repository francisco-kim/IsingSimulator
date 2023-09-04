using IsingMonteCarlo.Models;
using IsingMonteCarlo.Services;

var dimension = 2;
var latticeLength = 100;
var j = -1.0;
var h = 0.0;
//var spinUpdateMethod = SpinUpdateMethod.Glauber;
var spinUpdateMethod = SpinUpdateMethod.Wolff;
int? randomSeed = 41;

var boltzmannTemperature = 2.28;
string? filename = $"{latticeLength}_{boltzmannTemperature:0.0000}_1000000000.dat";
//string? filename = null;

var singleRunSimulation = new IsingSimulationSingleRun(filename,
                                                       dimension,
                                                       latticeLength,
                                                       boltzmannTemperature,
                                                       j,
                                                       h,
                                                       spinUpdateMethod,
                                                       randomSeed: randomSeed);

var thermalisationStepsInLatticeSizeUnit = 0;
if (filename is null)
{
    thermalisationStepsInLatticeSizeUnit = 100_000;
}

var iterationStepsBetweenMeasurements = singleRunSimulation.Simulation.TotalSpinsCount * 100;
const int measurementsCount = 50;
const int measurementsRepetitionCount = 10;

singleRunSimulation.RunWithMeasurements(iterationStepsBetweenMeasurements,
                                 measurementsCount,
                                 measurementsRepetitionCount,
                                 thermalisationStepsInLatticeSizeUnit,
                                 saveLattice: true,
                                 saveMeasurements: true);

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
