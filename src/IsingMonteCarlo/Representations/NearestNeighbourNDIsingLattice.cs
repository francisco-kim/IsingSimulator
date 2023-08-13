namespace IsingMonteCarlo.Representations;

public sealed class NearestNeighbourNDIsingLattice
{
    private const int LatticeSizeLowerBound = 3;

    private readonly List<int> _dNaryDigits;
    private readonly List<int> _decimalValuePerDNaryDigit;

    internal NearestNeighbourNDIsingLattice(
        int dimension,
        int latticeLength)
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


        Dimension = dimension;
        LatticeLength = latticeLength;

        _dNaryDigits = Enumerable.Range(0, Dimension)
                                 .Select(d => Convert.ToInt32(Math.Pow(LatticeLength, Dimension - d - 1)))
                                 .ToList();
        _decimalValuePerDNaryDigit = Enumerable.Range(0, Dimension)
                                               .Select(d => Convert.ToInt32(Math.Pow(LatticeLength, d)))
                                               .ToList();

        SpatialVectors = Enumerable.Range(0, TotalSpinsCount).Select(SpinIndexToSpatialVector).ToList();
        NeighboursIndices = Enumerable.Range(0, TotalSpinsCount).Select(FindNeighboursOfSpin).ToList();
    }

    public int Dimension { get; }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public List<List<int>> SpatialVectors { get; }

    public List<List<int>> NeighboursIndices { get; }

    // D-nary index: origin is on one of the corners of the hypercube
    public List<int> SpinIndexToSpatialPosition(int spinIndex)
    {
        if (spinIndex < 0 || spinIndex >= TotalSpinsCount)
        {
            throw new ArgumentOutOfRangeException(nameof(spinIndex),
                                                  "The spin index is not in the correct range.");
        }

        return SpinIndexToSpatialVector(spinIndex);
    }

    public int SpatialPositionToSpinIndex(IEnumerable<int> spatialPositionVector)
    {
        var spatialPosition = (spatialPositionVector?.ToList())
                              ?? throw new ArgumentNullException(nameof(spatialPositionVector));
        if (spatialPosition.Count > Dimension)
        {
            throw new ArgumentOutOfRangeException(nameof(spatialPositionVector),
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
