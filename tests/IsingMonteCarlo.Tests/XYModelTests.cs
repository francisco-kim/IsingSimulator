using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Tests;

public class XYModelTests
{
    [Fact]
    public void ColdStart_HasGroundStateEnergyMinusTwoPerSite()
    {
        var simulator = new XYModelSimulator(latticeLength: 16, hotStart: false);

        Assert.Equal(-2.0, simulator.GetEnergyPerSite(), precision: 12);
        Assert.Equal(1.0, simulator.GetMagnetisationMagnitude(), precision: 12);

        var (vortices, antivortices) = simulator.LocateVortices();
        Assert.Equal(0, vortices);
        Assert.Equal(0, antivortices);
    }

    [Fact]
    public void LowTemperature_StaysNearGroundStateWithSpinWaveCorrection()
    {
        var simulator = new XYModelSimulator(latticeLength: 32, hotStart: false, randomSeed: 41);
        var sweep = simulator.TotalSpinsCount;

        simulator.RunSteps(500L * sweep, beta: 1.0 / 0.2);

        // Equipartition of spin waves gives E/N ≈ −2 + T/2 = −1.9 at T = 0.2.
        var energyPerSite = simulator.GetEnergyPerSite();
        Assert.InRange(energyPerSite, -2.0, -1.8);
    }

    [Fact]
    public void Torus_VortexAndAntivortexCountsBalanceExactly()
    {
        var simulator = new XYModelSimulator(latticeLength: 32, hotStart: true, randomSeed: 7);
        var sweep = simulator.TotalSpinsCount;

        simulator.RunSteps(100L * sweep, beta: 1.0 / 1.4);

        var (vortices, antivortices) = simulator.LocateVortices();
        Assert.Equal(vortices, antivortices);
        Assert.True(vortices > 0, $"Expected free vortices at T = 1.4, found {vortices}.");
    }

    [Fact]
    public void KosterlitzThouless_VorticesProliferateAboveTheTransition()
    {
        var coldCount = TotalVortexCount(temperature: 0.4, seed: 11);
        var hotCount = TotalVortexCount(temperature: 1.5, seed: 12);

        // Below T_KT ≈ 0.893 vortices exist only as rare bound pairs; above,
        // they proliferate. Orders of magnitude apart — generous margin.
        Assert.True(
            hotCount > 10 * coldCount + 50,
            $"Vortex counts: T=0.4 -> {coldCount}, T=1.5 -> {hotCount}.");

        static int TotalVortexCount(double temperature, int seed)
        {
            var simulator = new XYModelSimulator(latticeLength: 32, hotStart: true, randomSeed: seed);
            var sweep = simulator.TotalSpinsCount;
            simulator.RunSteps(200L * sweep, beta: 1.0 / temperature);

            var total = 0;
            for (var snapshot = 0; snapshot < 10; snapshot++)
            {
                simulator.RunSteps(10L * sweep, beta: 1.0 / temperature);
                var (vortices, _) = simulator.LocateVortices();
                total += vortices;
            }

            return total;
        }
    }
}
