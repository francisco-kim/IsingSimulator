using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Services;

public static class Renormaliser
{
    private static readonly int RenormalisationSize = 3;
    private const int TwoPowNine = 512;
    private const int ThreePowSix = 729;

    public static void GenerateRenormalisedLatticeSpinConfigurationImages(
        NearestNeighbourNDIsingLattice<int> initialLattice,
        int finalLatticeSizeLimit,
        double? temperature = null,
        bool resize = true) =>
        GenerateRenormalisedLatticeSpinConfigurationImages(
            GenerateRenormalisedLattices(initialLattice, finalLatticeSizeLimit),
            temperature,
            resize);

    public static void GenerateRenormalisedLatticeSpinConfigurationImages(
        IEnumerable<NearestNeighbourNDIsingLattice<int>> lattices,
        double? temperature,
        bool resize = true)
    {
        var enumeratedLattices = lattices?.ToList() ?? throw new ArgumentNullException(nameof(lattices));
        var initialSpinConfigurations = enumeratedLattices.Select(l => l.Spins).ToList();
        var latticeLengths = enumeratedLattices.Select(l => l.LatticeLength).ToList();

        var bitmaps = initialSpinConfigurations.Select(DrawHelpers.GenerateGrayBitmapFrom2DList);

        var i = 0;
        foreach (var bitmap in bitmaps)
        {
            DrawHelpers.SaveBitmapAsPNG(
                bitmap,
                temperature is not null ? $"{latticeLengths[i]}_{temperature:0.00000}" : $"{latticeLengths[i]}",
                resize);
            ++i;
        }
    }

    public static List<NearestNeighbourNDIsingLattice<int>> GenerateRenormalisedLattices(
        NearestNeighbourNDIsingLattice<int> initialLattice,
        int finalLatticeSizeLimit)
    {
        var renormalisedLatticeConfigurations = new List<NearestNeighbourNDIsingLattice<int>> { initialLattice };
        while (initialLattice.LatticeLength >= finalLatticeSizeLimit)
        {
            initialLattice = PerformOneRenormalisationStepWithMajorityRule(initialLattice);
            renormalisedLatticeConfigurations.Add(initialLattice);
        }

        return renormalisedLatticeConfigurations;
    }

    public static List<List<int>> GenerateRenormalisedLatticeSpinConfigurationSpins(
        NearestNeighbourNDIsingLattice<int> initialLattice,
        int finalLatticeSizeLimit)
    {
        var renormalisedLatticeConfigurations = new List<List<int>>();
        while (initialLattice.LatticeLength >= finalLatticeSizeLimit)
        {
            initialLattice = PerformOneRenormalisationStepWithMajorityRule(initialLattice);
            renormalisedLatticeConfigurations.Add(initialLattice.Spins);
        }

        return renormalisedLatticeConfigurations;
    }

    private static NearestNeighbourNDIsingLattice<int> PerformOneRenormalisationStepWithMajorityRule(
        NearestNeighbourNDIsingLattice<int> initialLattice)
    {
        var initialLatticeLength = initialLattice.LatticeLength;
        if (initialLatticeLength % 3 is not 0)
        {
            throw new ArgumentException(nameof(initialLattice));
        }

        var initialLatticeConfiguration = initialLattice.Spins;

        var finalLatticeLength = initialLatticeLength / 3;
        var finalLatticeSize = finalLatticeLength * finalLatticeLength;
        var finalLatticeConfiguration = new List<int>(finalLatticeSize);

        for (var i = 0; i < finalLatticeSize; ++i)
        {
            var majorityRuleSpinSum = 0;
            var kernelReferenceCoordinate = initialLattice.SpinIndexToSpatialPosition(i * RenormalisationSize);
            for (var x = 0; x < RenormalisationSize; ++x)
            {
                for (var y = 0; y < RenormalisationSize; ++y)
                {
                    var spinIndexInKernel = initialLattice.SpatialPositionToSpinIndex(
                        new[] { x + kernelReferenceCoordinate[0], y + kernelReferenceCoordinate[1] });
                    majorityRuleSpinSum += initialLatticeConfiguration[spinIndexInKernel];
                }
            }

            var majorityRuleSpin = majorityRuleSpinSum > 0 ? 1 : -1;
            finalLatticeConfiguration.Add(majorityRuleSpinSum > 0 ? 1 : -1);
        }

        return new NearestNeighbourNDIsingLattice<int>(dimension: 2, finalLatticeLength, finalLatticeConfiguration);
    }
}
