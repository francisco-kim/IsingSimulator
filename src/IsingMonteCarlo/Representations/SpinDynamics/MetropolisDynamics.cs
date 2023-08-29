using IsingMonteCarlo.Services;

namespace IsingMonteCarlo.Representations.SpinDynamics;

public class MetropolisDynamics : ISpinDynamics
{
    private readonly IHamiltonian<int> _hamiltonian;
    private readonly int _totalSpinsCount;
    private readonly Random _random;

    public MetropolisDynamics(IHamiltonian<int> hamiltonian, int? randomSeed = null)
    {
        _hamiltonian = hamiltonian ?? throw new ArgumentNullException(nameof(hamiltonian));
        _totalSpinsCount = hamiltonian.Lattice.TotalSpinsCount;
        _random = randomSeed is not null ? new Random((int)randomSeed) : new Random();
    }

    public void FlipSpin(
        double beta,
        double j,
        double h,
        double? jY)
    {
        var chosenSite = _random.Next(minValue: 0, _totalSpinsCount);
        var flip = FlipWithMetropolis(randomProbability: _random.NextDouble(), chosenSite, beta, j, h, jY);

        if (flip)
        {
            _hamiltonian.FlipSpin(chosenSite);
        }
    }

    private bool FlipWithMetropolis(
        double randomProbability,
        int siteIndex,
        double beta,
        double j,
        double h,
        double? jY) =>
        randomProbability <= Math.Min(1.0, Math.Exp(-beta * _hamiltonian.GetDeltaEnergyOfSite(siteIndex, j, h, jY)));
}
