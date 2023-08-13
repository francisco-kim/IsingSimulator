using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Services;

var isSingleRun = false;

var dimension = 2;
var latticeLength = 100;

var monteCarlo = new MonteCarloSimulation(dimension, latticeLength);

// var boltzmannTemperature = 2.269;
var boltzmannTemperature = 2.5;
var beta = 1.0 / boltzmannTemperature;
var j = -1.0;
var h = 0.0;
var iterationLimit = monteCarlo.TotalSpinsCount * 1000;
var spinUpdateMethod = IsingMonteCarlo.Models.SpinUpdateMethod.Glauber;
var initialSpinDownRatio = 0.25;
var initialSpinConfiguration = SpinConfigurationBuilder.InitialiseLattice(monteCarlo.TotalSpinsCount,
                                                                          initialSpinDownRatio,
                                                                          randomSeed: 17);

if (isSingleRun)
{
    monteCarlo.RunMonteCarlo(beta,
                             j,
                             h,
                             iterationLimit,
                             spinUpdateMethod,
                             initialSpinConfiguration,
                             randomSeed: 17);

    var averageMagnetisation = monteCarlo.GetAverageMagnetisation();

    Console.WriteLine($"The average magnetisation: {averageMagnetisation}\n");
}
else
{
    var criticalRangeCount = 8;

    var temperaturesList = new List<double>(3) {0.0001, 0.8, 1.6};
    var criticalTemperatureRange = Enumerable.Range(0, criticalRangeCount)
                                             .Select(i => 0.5 / criticalRangeCount * i + 2.0);
    var highTemperatureRange = new List<double>(5) {2.5, 3.0, 3.5, 4.0, 4.5};

    temperaturesList.AddRange(criticalTemperatureRange);
    temperaturesList.AddRange(highTemperatureRange);

    var magnetisations = temperaturesList.Select((t, index)
                                                 => getMagnetisationFromEachRun(1.0 / t, index));

    var magnetisationData = temperaturesList.Zip(magnetisations);
    var magnetisationDataList = magnetisationData.Select(data => $"{{{data.First}, {data.Second}}}");
    Console.WriteLine(string.Join(", ", magnetisationDataList));

    double getMagnetisationFromEachRun(double beta, int iterationIndex)
    {
        monteCarlo.RunMonteCarlo(beta,
                                 j,
                                 h,
                                 iterationLimit,
                                 spinUpdateMethod,
                                 initialSpinConfiguration,
                                 randomSeed: 17);

        Console.WriteLine($"Iteration {iterationIndex} completed.");

        return monteCarlo.GetAverageMagnetisation();
    }
}
