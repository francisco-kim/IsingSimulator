using IsingMonteCarlo.Models;
using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Tests;

/// <summary>
///     Physics sanity checks with the repository's sign convention:
///     E_site = J * sum_neighbours(s_i * s_j) - h * s_i, so the ferromagnet has J = -1.
/// </summary>
public class PhysicsSanityTests
{
    private const double FerromagneticJ = -1.0;

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void OrderedLattice_HasGroundStateEnergyMinusDPerSite(int dimension)
    {
        var latticeLength = dimension is 1 ? 64 : 8;
        var simulator = new IsingMonteCarloSimulator(dimension, latticeLength);

        var energyPerSite = simulator.Hamiltonian.GetAverageEnergy(FerromagneticJ, h: 0.0);

        Assert.Equal(-dimension, energyPerSite, precision: 12);
    }

    [Fact]
    public void LowTemperature_ThermalisesToOrderedState()
    {
        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength: 16);
        var totalSpins = simulator.TotalSpinsCount;

        simulator.RunMonteCarlo(
            beta: 1.0 / 1.0,
            FerromagneticJ,
            h: 0.0,
            iterationLimit: 200L * totalSpins,
            SpinUpdateMethod.Wolff,
            randomSeed: 41);

        var absoluteMagnetisation = Math.Abs(simulator.Hamiltonian.GetAverageMagnetisation(FerromagneticJ, h: 0.0));

        // Exact Onsager magnetisation at T = 1 is > 0.999.
        Assert.True(
            absoluteMagnetisation > 0.95,
            $"|M| = {absoluteMagnetisation} after thermalisation at T = 1, expected near 1.");
    }

    [Fact]
    public void HighTemperature_MagnetisationNearZero()
    {
        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength: 32);
        var totalSpins = simulator.TotalSpinsCount;

        simulator.RunMonteCarlo(
            beta: 1.0 / 100.0,
            FerromagneticJ,
            h: 0.0,
            iterationLimit: 100L * totalSpins,
            SpinUpdateMethod.Glauber,
            randomSeed: 41);

        var absoluteMagnetisation = Math.Abs(simulator.Hamiltonian.GetAverageMagnetisation(FerromagneticJ, h: 0.0));

        Assert.True(
            absoluteMagnetisation < 0.1,
            $"|M| = {absoluteMagnetisation} at T = 100, expected near 0.");
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.7)]
    public void DeltaEnergy_MatchesRecomputedTotalEnergyDifference(double h)
    {
        var random = new Random(Seed: 17);
        var latticeLength = 8;
        var spins = Enumerable.Range(0, latticeLength * latticeLength)
                              .Select(_ => random.Next(2) is 0 ? 1 : -1)
                              .ToList();

        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength, spins);
        var hamiltonian = simulator.Hamiltonian;

        for (var trial = 0; trial < 20; trial++)
        {
            var site = random.Next(simulator.TotalSpinsCount);

            var energyBefore = FreshTotalEnergy(simulator, h);
            var predictedDelta = hamiltonian.GetDeltaEnergyOfSite(site, FerromagneticJ, h);

            hamiltonian.FlipSpin(site);
            var energyAfter = FreshTotalEnergy(simulator, h);

            Assert.Equal(energyAfter - energyBefore, predictedDelta, precision: 9);
        }

        // A fresh Hamiltonian recomputes the total energy from scratch instead
        // of using the incremental updates under test.
        static double FreshTotalEnergy(IsingMonteCarloSimulator sim, double h) =>
            new IsingHamiltonian(sim.Lattice).GetTotalEnergy(FerromagneticJ, h);
    }

    [Fact]
    public void TotalEnergy_CountsFieldTermOncePerSite()
    {
        // All-up 2D lattice: bonds contribute -2*N (J = -1, four bonds per site,
        // each shared by two sites), the field contributes -h*N.
        var latticeLength = 8;
        var totalSpins = latticeLength * latticeLength;
        var h = 0.5;

        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength);
        var totalEnergy = simulator.Hamiltonian.GetTotalEnergy(FerromagneticJ, h);

        Assert.Equal(-2.0 * totalSpins - h * totalSpins, totalEnergy, precision: 9);
    }

    [Fact]
    public void IncrementalMagnetisation_MatchesSpinSum()
    {
        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength: 16);
        simulator.RunMonteCarlo(
            beta: 1.0 / 2.5,
            FerromagneticJ,
            h: 0.0,
            iterationLimit: 20L * simulator.TotalSpinsCount,
            SpinUpdateMethod.Metropolis,
            randomSeed: 7);

        Assert.Equal(simulator.Lattice.Spins.Sum(), simulator.Hamiltonian.TotalMagnetisation, precision: 9);
    }

    [Fact]
    public void ResetObservables_MakesMeasurementBlocksIndependent()
    {
        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength: 8);
        var sweep = simulator.TotalSpinsCount;

        simulator.RunMonteCarlo(
            beta: 1.0 / 3.5,
            FerromagneticJ,
            h: 0.0,
            iterationLimit: 100L * sweep,
            SpinUpdateMethod.Glauber,
            randomSeed: 3);

        RunBlock(simulator);
        simulator.ResetObservables();
        var secondBlock = RunBlock(simulator);

        // The running average returned by the block must equal the mean of the
        // block's own samples; without the reset it would still include the
        // first block's samples.
        Assert.Equal(secondBlock.MagnetisationList.Average(), secondBlock.Magnetisation, precision: 12);
        Assert.Equal(secondBlock.EnergyList.Average(), secondBlock.Energy, precision: 12);

        BasicObservables RunBlock(IsingMonteCarloSimulator sim) =>
            sim.RunMonteCarloWithObservablesComputation(
                beta: 1.0 / 3.5,
                FerromagneticJ,
                h: 0.0,
                iterationsNeededForSingleChiXiMeasurement: 2 * sweep,
                measurementsCountForChiXiExpectationValue: 10,
                SpinUpdateMethod.Glauber);
    }

    [Fact]
    public void SpecificHeatAndSusceptibility_ArePositiveNearCriticality()
    {
        var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength: 16);
        var sweep = simulator.TotalSpinsCount;
        var beta = 1.0 / 2.269;

        simulator.RunMonteCarlo(
            beta,
            FerromagneticJ,
            h: 0.0,
            iterationLimit: 200L * sweep,
            SpinUpdateMethod.Wolff,
            randomSeed: 5);

        var observables = simulator.RunMonteCarloWithObservablesComputation(
            beta,
            FerromagneticJ,
            h: 0.0,
            iterationsNeededForSingleChiXiMeasurement: 5 * sweep,
            measurementsCountForChiXiExpectationValue: 40,
            SpinUpdateMethod.Wolff);

        Assert.True(observables.SpecificHeat > 0.0, $"C = {observables.SpecificHeat}");
        Assert.True(observables.Susceptibility > 0.0, $"Chi = {observables.Susceptibility}");
    }

    [Fact]
    public void WolffAndGlauber_AgreeOnEnergyAtModerateTemperature()
    {
        // Cross-check of the two dynamics: both must sample the same equilibrium
        // distribution. T = 3.5 (paramagnetic phase, short correlation time).
        var beta = 1.0 / 3.5;

        var wolffEnergy = EquilibriumEnergy(SpinUpdateMethod.Wolff, beta, seed: 11);
        var glauberEnergy = EquilibriumEnergy(SpinUpdateMethod.Glauber, beta, seed: 12);

        Assert.True(
            Math.Abs(wolffEnergy - glauberEnergy) < 0.05,
            $"E_Wolff = {wolffEnergy}, E_Glauber = {glauberEnergy}");

        static double EquilibriumEnergy(SpinUpdateMethod method, double beta, int seed)
        {
            var simulator = new IsingMonteCarloSimulator(dimension: 2, latticeLength: 16);
            var sweep = simulator.TotalSpinsCount;

            simulator.RunMonteCarlo(
                beta,
                FerromagneticJ,
                h: 0.0,
                iterationLimit: 200L * sweep,
                method,
                randomSeed: seed);

            var observables = simulator.RunMonteCarloWithObservablesComputation(
                beta,
                FerromagneticJ,
                h: 0.0,
                iterationsNeededForSingleChiXiMeasurement: 2 * sweep,
                measurementsCountForChiXiExpectationValue: 100,
                method);

            return observables.Energy;
        }
    }
}
