using System.Numerics;

namespace IsingMonteCarlo.Representations;

public sealed class NearestNeighbourNDIsingLattice<T> where T : INumber<T>
{
    private const int LatticeSizeLowerBound = 3;

    private readonly List<int> _dNaryDigits;
    private readonly List<int> _decimalValuePerDNaryDigit;

    internal NearestNeighbourNDIsingLattice(
        int dimension,
        int latticeLength,
        IEnumerable<T>? initialSpinConfiguration = null)
    {
        if (dimension < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(dimension),
                $"The dimension must be greater than 1, but {dimension} was given.");
        }

        if (latticeLength < LatticeSizeLowerBound)
        {
            throw new ArgumentOutOfRangeException(
                nameof(latticeLength),
                $"The lattice length be greater than {LatticeSizeLowerBound}, but {latticeLength} was given.");
        }

        TotalSpinsCount = Convert.ToInt32(Math.Pow(latticeLength, dimension));

        var initialSpins = initialSpinConfiguration is not null
                               ? initialSpinConfiguration.ToList()
                               : Enumerable.Repeat(T.MultiplicativeIdentity, TotalSpinsCount).ToList();

        if (initialSpins.Count != TotalSpinsCount)
        {
            throw new ArgumentException(
                $"There must be (lattice length)^(dimension) spins, but {initialSpins.Count} spins were given.",
                nameof(initialSpinConfiguration));
        }

        if (initialSpins.Any(
                spin => spin != T.MultiplicativeIdentity
                     && spin != -T.MultiplicativeIdentity))
        {
            throw new ArgumentException(
                "The spins must have integer values +1 or -1.",
                nameof(initialSpinConfiguration));
        }

        Dimension = dimension;
        LatticeLength = latticeLength;

        _dNaryDigits = Enumerable.Range(start: 0, Dimension)
                                 .Select(d => Convert.ToInt32(Math.Pow(LatticeLength, Dimension - d - 1)))
                                 .ToList();
        _decimalValuePerDNaryDigit = Enumerable.Range(start: 0, Dimension)
                                               .Select(d => Convert.ToInt32(Math.Pow(LatticeLength, d)))
                                               .ToList();

        SpatialVectors = Enumerable.Range(start: 0, TotalSpinsCount).Select(SpinIndexToSpatialVector).ToList();
        NeighboursIndices = Enumerable.Range(start: 0, TotalSpinsCount).Select(FindNeighboursOfSpin).ToList();
        Spins = initialSpins;
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public List<T> Spins { get; set; }

    public List<List<int>> SpatialVectors { get; }

    public List<List<int>> NeighboursIndices { get; }

    // D-nary index: origin is on one of the corners of the hypercube
    public List<int> SpinIndexToSpatialPosition(int spinIndex)
    {
        if (spinIndex < 0 || spinIndex >= TotalSpinsCount)
        {
            throw new ArgumentOutOfRangeException(
                nameof(spinIndex),
                "The spin index is not in the correct range.");
        }

        return SpinIndexToSpatialVector(spinIndex);
    }

    public int SpatialPositionToSpinIndex(IEnumerable<int> spatialPositionVector)
    {
        var spatialPosition = spatialPositionVector?.ToList()
                           ?? throw new ArgumentNullException(nameof(spatialPositionVector));
        if (spatialPosition.Count > Dimension)
        {
            throw new ArgumentOutOfRangeException(
                nameof(spatialPositionVector),
                "The spatial position vector must have the dimension of the Ising model.");
        }

        return SpatialVectorToSpinIndex(spatialPosition);
    }

    // D-nary index: origin is on one of the corners of the hypercube
    internal List<int> SpinIndexToSpatialVector(int spinIndex)
    {
        var listOfPositions = new List<int>(Dimension);
        var remainder = spinIndex;
        for (var dimensionIndex = 0; dimensionIndex < Dimension; dimensionIndex++)
        {
            var (quotient, newRemainder) = Math.DivRem(remainder, _dNaryDigits[dimensionIndex]);
            listOfPositions.Add(quotient);
            remainder = newRemainder;
        }

        listOfPositions.Reverse();

        return listOfPositions;
    }

    // D-nary-to-decimal conversion
    internal int SpatialVectorToSpinIndex(List<int> dNaryNumber) =>
        dNaryNumber.Select((digit, index) => digit * _decimalValuePerDNaryDigit[index])
                   .Sum();

    private List<int> FindNeighboursOfSpin(int spinIndex)
    {
        var spatialVector = SpatialVectors[spinIndex];

        var listOfNeighbours = new List<int>(2 * Dimension);
        for (var dimensionIndex = 0; dimensionIndex < Dimension; dimensionIndex++)
        {
            var positiveDirectionNeighbour = spatialVector.ToList();
            positiveDirectionNeighbour[dimensionIndex] = (spatialVector[dimensionIndex] + 1) % LatticeLength;

            // Do this instead of Modulo(index - 1, LatticeLength) to avoid minus-sign problem
            var negativeDirectionNeighbour = spatialVector.ToList();
            negativeDirectionNeighbour[dimensionIndex] = (spatialVector[dimensionIndex] + (LatticeLength - 1))
                                                       % LatticeLength;

            listOfNeighbours.Add(SpatialVectorToSpinIndex(positiveDirectionNeighbour));
            listOfNeighbours.Add(SpatialVectorToSpinIndex(negativeDirectionNeighbour));
        }

        return listOfNeighbours;
    }
}
