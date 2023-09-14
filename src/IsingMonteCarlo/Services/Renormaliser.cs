using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Services;

public static class Renormaliser
{
    private static readonly int RenormalisationSize = 3;

    public static void GenerateRenormalisedLatticeSpinConfigurationImages(
        IEnumerable<int> initialSpinConfiguration,
        int finalLatticeSizeLimit,
        double temperature = double.NaN,
        bool resize = true)
    {
        var startingLattice = initialSpinConfiguration.ToList()
                           ?? throw new ArgumentNullException(nameof(initialSpinConfiguration));

        var numberedFolderName = GetSubFolderNumberedName("renormalisation").ToArray();
        GenerateLatticeSpinConfigurationImage(
            startingLattice,
            temperature,
            numberedFolderName,
            resize);

        var latticeLength = Convert.ToInt32(Math.Sqrt(startingLattice.Count));
        while (latticeLength >= finalLatticeSizeLimit)
        {
            startingLattice = GenerateRenormalisedLattice(startingLattice, finalLatticeSizeLimit);
            latticeLength = Convert.ToInt32(Math.Sqrt(startingLattice.Count));
            GenerateLatticeSpinConfigurationImage(
                startingLattice,
                temperature,
                numberedFolderName,
                resize);
        }
    }

    public static void GenerateZoomedLatticeSpinConfigurationImages(
        IEnumerable<int> initialSpinConfiguration,
        int finalLatticeSizeLimit,
        double temperature = double.NaN,
        bool resize = true)
    {
        var startingLattice = initialSpinConfiguration.ToList()
                           ?? throw new ArgumentNullException(nameof(initialSpinConfiguration));

        var numberedFolderName = GetSubFolderNumberedName("zoom").ToArray();
        GenerateLatticeSpinConfigurationImage(
            startingLattice,
            temperature,
            numberedFolderName,
            resize);

        var latticeLength = Convert.ToInt32(Math.Sqrt(startingLattice.Count));
        while (latticeLength >= finalLatticeSizeLimit && latticeLength % 3 == 0)
        {
            startingLattice = GenerateZoomedLattice(startingLattice, finalLatticeSizeLimit);
            latticeLength = Convert.ToInt32(Math.Sqrt(startingLattice.Count));
            GenerateLatticeSpinConfigurationImage(
                startingLattice,
                temperature,
                numberedFolderName,
                resize);
        }
    }

    public static void GenerateLatticeSpinConfigurationImage(
        IEnumerable<int> initialSpinConfiguration,
        double temperature,
        IEnumerable<string> folderName,
        bool resize = true)
    {
        var latticeSpinConfiguration = initialSpinConfiguration.ToList()
                                    ?? throw new ArgumentNullException(nameof(initialSpinConfiguration));
        var latticeLength = Convert.ToInt32(Math.Sqrt(latticeSpinConfiguration.Count));

        var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(latticeSpinConfiguration);

        DrawHelpers.SaveBitmapAsPNGInSpecifiedFolder(
            bitmap,
            latticeLength,
            temperature,
            folderName,
            resize: resize);
    }

    public static List<int> GenerateRenormalisedLattice(
        IEnumerable<int> initialSpinConfiguration,
        int finalLatticeSizeLimit) =>
        PerformOneRenormalisationStepWithMajorityRule(
            Fast2DIsingMonteCarloSimulator.Get2DArrayFrom1DArray(initialSpinConfiguration));

    public static List<int> GenerateZoomedLattice(
        IEnumerable<int> initialSpinConfiguration,
        int finalLatticeSizeLimit) =>
        PerformZoom(Fast2DIsingMonteCarloSimulator.Get2DArrayFrom1DArray(initialSpinConfiguration));

    private static List<int> PerformOneRenormalisationStepWithMajorityRule(
        int[][] initialLattice)
    {
        var initialLatticeLength = initialLattice.Length;

        if (initialLatticeLength % 3 is not 0)
        {
            throw new ArgumentException(nameof(initialLattice));
        }

        //initialLattice = null;
        //GC.Collect();
        //GC.WaitForPendingFinalizers();

        var finalLatticeLength = initialLatticeLength / 3;
        var finalLatticeSize = finalLatticeLength * finalLatticeLength;
        var finalLatticeConfiguration = new List<int>(finalLatticeSize);

        for (var i = 0; i < finalLatticeSize; i++)
        {
            var majorityRuleSpinSum = 0;
            var (quotient, kernelReferenceCoordinateX) = Math.DivRem(i * RenormalisationSize, initialLatticeLength);
            var kernelReferenceCoordinateY = quotient * RenormalisationSize;
            for (var x = 0; x < RenormalisationSize; x++)
            {
                for (var y = 0; y < RenormalisationSize; y++)
                {
                    majorityRuleSpinSum +=
                        initialLattice[x + kernelReferenceCoordinateX][y + kernelReferenceCoordinateY];
                }
            }

            var majorityRuleSpin = majorityRuleSpinSum > 0 ? 1 : -1;
            finalLatticeConfiguration.Add(majorityRuleSpin);
        }

        return finalLatticeConfiguration;
    }

    private static List<int> PerformZoom(
        int[][] initialLattice)
    {
        var initialLatticeLength = initialLattice.Length;

        //initialLattice = null;
        //GC.Collect();
        //GC.WaitForPendingFinalizers();

        var finalLatticeLength = initialLatticeLength / 3 * 2;
        var finalLatticeSize = finalLatticeLength * finalLatticeLength;
        var finalLatticeConfiguration = new List<int>(finalLatticeSize);

        for (var y = 0; y < finalLatticeLength; y++)
        {
            for (var x = 0; x < finalLatticeLength; x++)
            {
                finalLatticeConfiguration.Add(initialLattice[x][y]);
            }
        }

        return finalLatticeConfiguration;
    }

    private static IEnumerable<string> GetSubFolderNumberedName(string folderName)
    {
        var newFolderName = new string[] { };
        for (var folderIndex = 0; folderIndex < 999; folderIndex++)
        {
            var imagesDirectory =
                FileHelpers.GetDataRootDirectory(new[] { folderName, folderIndex.ToString(format: "D3") });

            if (Directory.Exists(imagesDirectory))
            {
                continue;
            }

            newFolderName = new[] { folderName, folderIndex.ToString(format: "D3") };
            break;
        }

        return newFolderName;
    }
}
