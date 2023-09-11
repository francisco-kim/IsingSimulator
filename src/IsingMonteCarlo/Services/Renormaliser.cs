using IsingMonteCarlo.Helpers;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Services;

public static class Renormaliser
{
    private static readonly int RenormalisationSize = 3;

    public static void GenerateRenormalisedLatticeSpinConfigurationImages(
        NearestNeighbourNDIsingLattice<int> initialLattice,
        int finalLatticeSizeLimit,
        double temperature = double.NaN,
        bool resize = true)
    {
        var startingLattice = initialLattice;
        while (startingLattice.LatticeLength >= finalLatticeSizeLimit)
        {
            startingLattice = GenerateRenormalisedLattice(initialLattice, finalLatticeSizeLimit);
            GenerateRenormalisedLatticeSpinConfigurationImage(
                startingLattice,
                temperature,
                resize);
        }
    }

    public static void GenerateRenormalisedLatticeSpinConfigurationImage(
        NearestNeighbourNDIsingLattice<int> lattice,
        double temperature,
        bool resize = true)
    {
        var initialSpinConfiguration = lattice.Spins;
        var latticeLength = lattice.LatticeLength;

        var bitmap = DrawHelpers.GenerateGrayBitmapFrom2DList(initialSpinConfiguration);

        DrawHelpers.SaveBitmapAsPNG(
            bitmap,
            latticeLength,
            temperature,
            resize);
    }

    public static NearestNeighbourNDIsingLattice<int> GenerateRenormalisedLattice(
        NearestNeighbourNDIsingLattice<int> initialLattice,
        int finalLatticeSizeLimit) =>
        PerformOneRenormalisationStepWithMajorityRule(initialLattice);

    public static List<int> GenerateRenormalisedLatticeSpinConfigurationSpins(
        NearestNeighbourNDIsingLattice<int> initialLattice,
        int finalLatticeSizeLimit) =>
        PerformOneRenormalisationStepWithMajorityRule(initialLattice).Spins;

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
