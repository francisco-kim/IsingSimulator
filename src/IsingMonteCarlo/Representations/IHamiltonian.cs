using System.Numerics;

namespace IsingMonteCarlo.Representations;

public interface IHamiltonian<T> where T : INumber<T>
{
    NearestNeighbourNDIsingLattice<T> Lattice { get; }

    double TotalEnergy { get; }

    double TotalMagnetisation { get; }

    double GetTotalEnergy(
        double j,
        double h,
        double? jY = null);

    void FlipSpin(int spinIndex);

    public void FlipSpinWithPropertiesUpdate(
        int spinIndex,
        double j,
        double h,
        double? jY = null);

    /// <summary>
    ///     Calculates the average energy per site.
    ///     Its absolute value is the hypercube dimension D
    ///     (with an isotropic J-coupling and without an external field).
    /// </summary>
    double GetAverageEnergy(
        double j,
        double h,
        double? jY = null);

    double GetAverageMagnetisation(
        double j,
        double h,
        double? jY = null);

    double GetDeltaEnergyOfSite(
        int spinIndex,
        double j,
        double h,
        double? jY = null);

    double GetEnergyOfSite(
        int spinIndex,
        double j,
        double h,
        double? jY = null);
}
