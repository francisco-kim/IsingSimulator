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
        TotalMagnetisation = Lattice.Spins.Sum();
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
            // Each bond is counted twice when summing per-site bond energies,
            // but the field term -h*s_i is counted once per site.
            var bondEnergySum = 0.0;
            for (var spinIndex = 0; spinIndex < _totalSpinsCount; spinIndex++)
            {
                bondEnergySum += GetBondEnergyOfSite(spinIndex, j, jY);
            }

            TotalEnergy = 0.5 * bondEnergySum - h * TotalMagnetisation;
        }

        return TotalEnergy;
    }

    public void FlipSpin(int spinIndex) => Lattice.Spins[spinIndex] *= -1;

    public void InvalidateTotalEnergy() => TotalEnergy = double.NaN;

    public void FlipSpinWithPropertiesUpdate(
        int spinIndex,
        double j,
        double h,
        double? jY = null)
    {
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
        double? jY = null) => TotalMagnetisation / _totalSpinsCount;

    public double GetDeltaEnergyOfSite(
        int spinIndex,
        double j,
        double h,
        double? jY = null) => -2.0 * GetEnergyOfSite(spinIndex, j, h, jY);

    public double GetEnergyOfSite(
        int spinIndex,
        double j,
        double h,
        double? jY = null) =>
        GetBondEnergyOfSite(spinIndex, j, jY) - h * Lattice.Spins[spinIndex];

    private double GetBondEnergyOfSite(
        int spinIndex,
        double j,
        double? jY)
    {
        var spins = Lattice.Spins;
        var neighbours = _neighboursIndices[spinIndex];
        var spinValue = spins[spinIndex];

        if (jY is null)
        {
            var neighbourSum = 0;
            for (var i = 0; i < neighbours.Count; i++)
            {
                neighbourSum += spins[neighbours[i]];
            }

            return j * spinValue * neighbourSum;
        }

        if (_dimension != 2)
        {
            throw new ArgumentException(
                "The coupling constant $J_{Y}$ is valid only in the 2D Ising model.",
                nameof(jY));
        }

        // Neighbour ordering per site is [+x, -x, +y, -y].
        var xBondsSum = spins[neighbours[0]] + spins[neighbours[1]];
        var yBondsSum = spins[neighbours[2]] + spins[neighbours[3]];

        return spinValue * (j * xBondsSum + (double)jY * yBondsSum);
    }
}
