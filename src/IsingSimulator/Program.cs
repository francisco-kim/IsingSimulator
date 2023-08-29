using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;
using IsingMonteCarlo.Services;

var isSingleRun = true;

var dimension = 2;
var latticeLength = 100;

var totalSpinsCount = Convert.ToInt32(Math.Pow(latticeLength, dimension));

var initialSpinDownRatio = 0.25;
var initialSpinConfiguration = SpinConfigurationBuilder.InitialiseLattice(totalSpinsCount,
                                                                          initialSpinDownRatio,
                                                                          randomSeed: 17);
var monteCarlo = new IsingMonteCarloSimulation(dimension, latticeLength, initialSpinConfiguration);

// var boltzmannTemperature = 2.269;
var boltzmannTemperature = 1.8;
var beta = 1.0 / boltzmannTemperature;
var j = -1.0;
var h = 0.0;
var iterationLimit = monteCarlo.TotalSpinsCount * 1000;
var spinUpdateMethod = SpinUpdateMethod.Wolff;

if (isSingleRun)
{
    monteCarlo.RunMonteCarlo(beta,
                             j,
                             h,
                             iterationLimit,
                             spinUpdateMethod,
                             randomSeed: 17);

    var averageMagnetisation = monteCarlo.Hamiltonian.GetAverageMagnetisation(j, h);

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
                                 randomSeed: 17);

        Console.WriteLine($"Iteration {iterationIndex} completed.");

        return new IsingHamiltonian(monteCarlo.Lattice).GetAverageMagnetisation(j, h);
    }
}
