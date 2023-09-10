using IsingMonteCarlo.Models;

namespace IsingMonteCarlo.Helpers;

public static class SpinConfigurationBuilder
{
    public static InitialLatticeData InitialiseLattice(string? filename,
                                              int dimension = 2,
                                              int latticeLength = 100,
                                              double initialSpinDownRatio = 0.25,
                                              int? randomSeed = 41)
    {
        if (filename is not null)
        {
            var fullPathWithFilename = FileHelpers.GetFullPathWithFilename(filename);

            var spinConfiguration =  FileHelpers.LoadSpinConfiguration(fullPathWithFilename,
                                                         out var boltzmannTemperature,
                                                         out var previousIterationCountIn100MCSweepUnit);

            return (spinConfiguration, boltzmannTemperature, previousIterationCountIn100MCSweepUnit);
        }

        var totalSpinsCount = Convert.ToInt32(Math.Pow(latticeLength, dimension));

        return (InitialiseRandomLattice(totalSpinsCount,
                                       initialSpinDownRatio,
                                       randomSeed), boltzmannTemperature: double.NaN, previousIterationCount: 0);
    }

    public static List<int> InitialiseRandomLattice(int totalSpinsCount,
                                                    double initialSpinDownRatio,
                                                    int? randomSeed = null)
    {
        if (totalSpinsCount < 1)
        {
            throw new ArgumentException("There cannot be zero or less spins.", nameof(totalSpinsCount));
        }

        if (initialSpinDownRatio < 0.0 || initialSpinDownRatio > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(initialSpinDownRatio),
                $"The initial spin-down ratio must be between 0 and 1, but {initialSpinDownRatio} was given.");
        }

        var random = randomSeed is not null ? new Random((int)randomSeed) : new Random();

        var spinsToSet = new List<int>(totalSpinsCount);
        for (var i = 0; i < totalSpinsCount; i++)
        {
            spinsToSet.Add(random.NextDouble() > initialSpinDownRatio ? 1 : -1);
        }

        return spinsToSet;
    }
}
