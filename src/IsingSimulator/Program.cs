using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Models;
using IsingMonteCarlo.Services;

var loadFile = true;
loadFile = false;
var isSingleRun = true;

var dimension = 2;
var latticeLength = 100;

var totalSpinsCount = Convert.ToInt32(Math.Pow(latticeLength, dimension));

var previousIterationCount = 0;
var boltzmannTemperature = 0.0;
var initialSpinConfiguration = new List<int>(totalSpinsCount);
if (loadFile)
{
    var workingDirectory = Directory.GetCurrentDirectory();
    var rootDirectory = Directory.GetParent(workingDirectory)?.Parent ?? throw new ArgumentException(nameof(workingDirectory));
    var dataDirectory = Path.GetFullPath(Path.Combine(rootDirectory.FullName, "data"));
    var firstFileName = new DirectoryInfo(dataDirectory).EnumerateFiles()
                                                        .Select(file => file.FullName)
                                                        .FirstOrDefault();

    initialSpinConfiguration = LatticeConfigurationSaver.LoadLattice(firstFileName, out boltzmannTemperature, out previousIterationCount);
}
else
{
    var initialSpinDownRatio = 0.25;
    initialSpinConfiguration = SpinConfigurationBuilder.InitialiseLattice(totalSpinsCount,
                                                                              initialSpinDownRatio,
                                                                              randomSeed: 17);
    var temperature = 0.0;
    boltzmannTemperature = 4.0;
}

var monteCarlo = new IsingMonteCarloSimulation(dimension, latticeLength, initialSpinConfiguration);

// var boltzmannTemperature = 2.269;
var beta = 1.0 / boltzmannTemperature;
var j = -1.0;
var h = 0.0;
var thermalisationCount = monteCarlo.TotalSpinsCount * 20_000;
var iterationStepsBetweenMeasurements = monteCarlo.TotalSpinsCount * 50;
var measurementsCount = 100;
var spinUpdateMethod = SpinUpdateMethod.Wolff;

if (isSingleRun)
{
    // Thermalisation
    monteCarlo.RunMonteCarlo(
        beta,
        j,
        h,
        thermalisationCount,
        spinUpdateMethod,
        randomSeed: 17);

    var measuredMagnetisation = monteCarlo.Hamiltonian.GetAverageMagnetisation(j, h);
    Console.WriteLine($"The measured magnetisation: {measuredMagnetisation}\n\n");

    monteCarlo.RunMonteCarloWithObservablesComputation(
        beta,
        j,
        h,
        iterationStepsBetweenMeasurements,
        measurementsCount,
        spinUpdateMethod,
        randomSeed: 17);

    Console.WriteLine($" M  = {monteCarlo.Magnetisation} +- {monteCarlo.MagnetisationSigma}");
    Console.WriteLine($"M^2 = {monteCarlo.MagnetisationSquared} +- {monteCarlo.MagnetisationSquaredSigma}");
    Console.WriteLine($"|M| = {monteCarlo.MagnetisationAbsolute} +- {monteCarlo.MagnetisationAbsoluteSigma}");
    Console.WriteLine($"Chi = {monteCarlo.Susceptibility}");
    Console.WriteLine($" Xi = {monteCarlo.RenormalisedCorrelationLength}");
    // var correlationLengthListForPrintingX = correlationLengthList.Where(xi => xi.InXDirection is not double.NaN)
    //                                                              .Select(xi => $"{xi.InXDirection}");
    // var correlationLengthListForPrintingY = correlationLengthList.Where(xi => xi.InYDirection is not double.NaN)
    //                                                              .Select(xi => $"{xi.InYDirection}");
    // Console.WriteLine(
    //     string.Join(",", correlationLengthListForPrintingX.Concat(correlationLengthListForPrintingY).ToList()));

    LatticeConfigurationSaver.SaveLattice(
        monteCarlo.Lattice.Spins,
        boltzmannTemperature,
        previousIterationCount + thermalisationCount + measurementsCount * iterationStepsBetweenMeasurements,
        isBinary: true);
}
else
{
    var criticalRangeCount = 8;

    var temperaturesList = new List<double>(3) { 0.0001, 0.8, 1.6 };
    var criticalTemperatureRange = Enumerable.Range(0, criticalRangeCount)
                                             .Select(i => 0.5 / criticalRangeCount * i + 2.0);
    var highTemperatureRange = new List<double>(5) { 2.5, 3.0, 3.5, 4.0, 4.5 };

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
                                 thermalisationCount,
                                 spinUpdateMethod,
                                 randomSeed: 17);

        Console.WriteLine($"Iteration {iterationIndex} completed.");

        return monteCarlo.Hamiltonian.GetAverageMagnetisation(j, h);
    }
}
