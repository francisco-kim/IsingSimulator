namespace IsingMonteCarlo.Representations;

public sealed class IsingHamiltonian : IHamiltonian<int>
{
    private readonly List<List<int>> _neighboursIndices;
    private readonly int _dimension;
    private readonly int _totalSpinsCount;

    public IsingHamiltonian(NearestNeighbourNDIsingLattice<int> lattice)
    {
        if (lattice is null)
        {
            throw new ArgumentNullException(nameof(lattice));
        }

        Lattice = lattice;
        _dimension = lattice.Dimension;
        _neighboursIndices = lattice.NeighboursIndices;
        _totalSpinsCount = lattice.TotalSpinsCount;

        TotalEnergy = double.NaN;
        TotalMagnetisation = int.MinValue;
    }

    public NearestNeighbourNDIsingLattice<int> Lattice { get; }

    public double TotalEnergy { get; private set; }

    public double TotalMagnetisation { get; private set; }

    public double GetTotalEnergy(
        double j,
        double h,
        double? jY = null)
    {
        if (TotalEnergy is double.NaN)
        {
            TotalEnergy = 0.5 * Enumerable.Range(0, _totalSpinsCount)
                        .Select(i => GetEnergyOfSite(i, j, h, jY))
                        .Sum();
        }

        return TotalEnergy;
    }

    public void FlipSpin(int spinIndex) => Lattice.Spins[spinIndex] *= -1;

    public void FlipSpinWithPropertiesUpdate(
        int spinIndex,
        double j,
        double h,
        double? jY = null)
    {
        var totalEnergy = GetTotalEnergy(j, h, jY);
        var totalMagnetisation = GetTotalMagnetisation();

        TotalEnergy += GetDeltaEnergyOfSite(spinIndex, j, h, jY);
        TotalMagnetisation += -2 * Lattice.Spins[spinIndex];

        Lattice.Spins[spinIndex] *= -1;
    }

    /// <summary>
    ///     Calculates the average energy per site.
    ///     Its absolute value is the hypercube dimension D
    ///     (with an isotropic J-coupling and without an external field).
    /// </summary>
    public double GetAverageEnergy(
        double j,
        double h,
        double? jY = null) => GetTotalEnergy(j, h, jY) / _totalSpinsCount;

    public double GetAverageMagnetisation(
        double j,
        double h,
        double? jY = null) => (double)GetTotalMagnetisation() / _totalSpinsCount;

    public double GetDeltaEnergyOfSite(
        int spinIndex,
        double j,
        double h,
        double? jY = null) => -2.0 * GetEnergyOfSite(spinIndex, j, h, jY);

    public double GetEnergyOfSite(
        int spinIndex,
        double j,
        double h,
        double? jY = null)
    {
        if (jY is not null && _dimension != 2)
        {
            throw new ArgumentException(
                "The coupling constant $J_{Y}$ is valid only in the 2D Ising model.",
                nameof(jY));
        }
        var spinValue = Lattice.Spins[spinIndex];

        // Possible in 2D only
        if (jY is not null && _dimension is 2)
        {
            var xBonds = j * _neighboursIndices[spinIndex].Take(2)
                                                          .Select(
                                                            neighbourIndex => 
                                                                spinValue
                                                                * Lattice.Spins[neighbourIndex])
                                                          .Sum();
            var yBonds = jY * _neighboursIndices[spinIndex].Skip(2).Take(2)
                                                           .Select(
                                                            neighbourIndex => 
                                                                spinValue
                                                                * Lattice.Spins[neighbourIndex])
                                                            .Sum();
            
            return xBonds + (double)yBonds - h * spinValue;
        }

        return j * _neighboursIndices[spinIndex].Select(neighbourIndex => Lattice.Spins[neighbourIndex]
                                                                          * spinValue)
                                                .Sum()
               - h * spinValue;
    }

    private double GetTotalMagnetisation()
    {
        if (TotalMagnetisation is int.MinValue)
        {
            TotalMagnetisation = Lattice.Spins.Sum();
        }

        return TotalMagnetisation;
    }
}
