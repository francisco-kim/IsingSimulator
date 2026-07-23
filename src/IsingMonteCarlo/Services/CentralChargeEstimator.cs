namespace IsingMonteCarlo.Services;

/// <summary>
///     Exact extraction of the central charge and leading operator dimensions of the
///     2D square-lattice Ising critical point, via Kaufman's free-fermion
///     diagonalisation of the transfer matrix. No Monte Carlo is involved: every
///     quantity here is a closed-form sum or integral over the fermion dispersion,
///     evaluated on a finite cylinder of circumference L.
/// </summary>
public static class CentralChargeEstimator
{
    /// <summary>Critical coupling K_c = J / k_B T_c, fixed by sinh(2K_c) = 1.</summary>
    public static readonly double CriticalCoupling = 0.5 * Math.Log(1.0 + Math.Sqrt(2.0));

    /// <summary>
    ///     Free-fermion dispersion relation cosh γ(q) = cosh(2K)·coth(2K) − cos q.
    ///     At K = K_c this reduces to cosh γ(q) = 2 − cos q, and γ(q) → v|q| as
    ///     q → 0 with speed v = 1 for the isotropic lattice.
    /// </summary>
    public static double Gamma(double q, double coupling)
    {
        var cosh2K = Math.Cosh(2.0 * coupling);
        var sinh2K = Math.Sinh(2.0 * coupling);
        var argument = cosh2K * cosh2K / sinh2K - Math.Cos(q);
        return Math.Acosh(Math.Max(1.0, argument));
    }

    /// <summary>
    ///     Reduced free energy per site, g(L) = −f(L), of the width-L cylinder's
    ///     largest transfer-matrix eigenvalue (Neveu–Schwarz / antiperiodic sector,
    ///     the sector containing the ground state).
    /// </summary>
    public static double FreeEnergyG(int latticeWidth, double coupling)
    {
        var sinh2K = Math.Sinh(2.0 * coupling);
        var sum = 0.0;
        for (var r = 0; r < latticeWidth; r++)
        {
            var q = Math.PI * (2 * r + 1) / latticeWidth;
            sum += Gamma(q, coupling);
        }

        return 0.5 * Math.Log(2.0 * sinh2K) + sum / (2.0 * latticeWidth);
    }

    /// <summary>
    ///     Bulk (L → ∞) reduced free energy per site, by quadrature of the same sum
    ///     turned into an integral over the Brillouin zone.
    /// </summary>
    public static double FreeEnergyGInfinity(double coupling, int quadraturePoints = 20_000)
    {
        var sinh2K = Math.Sinh(2.0 * coupling);
        var sum = 0.0;
        for (var i = 0; i < quadraturePoints; i++)
        {
            var q = 2.0 * Math.PI * (i + 0.5) / quadraturePoints;
            sum += Gamma(q, coupling);
        }

        return 0.5 * Math.Log(2.0 * sinh2K) + 0.5 * (sum / quadraturePoints);
    }

    /// <summary>
    ///     Effective central charge from the finite-size free-energy amplitude at
    ///     coupling K: c_eff(L) = 6 L² (g(L) − g_∞) / π. Equals the true central
    ///     charge only at the critical point; decays towards 0 away from it as the
    ///     fermion gap γ(0) = 2|K − K_c| opens.
    /// </summary>
    public static double EffectiveCentralCharge(int latticeWidth, double coupling)
    {
        var g = FreeEnergyG(latticeWidth, coupling);
        var gInfinity = FreeEnergyGInfinity(coupling);
        return 6.0 * latticeWidth * latticeWidth * (g - gInfinity) / Math.PI;
    }

    /// <summary>
    ///     Least-squares fit at K_c of g(L) = g_∞ + A/L² + B/L⁴ over the given
    ///     widths. The central charge is c = 6A/π (Cardy's formula, speed v = 1 for
    ///     the isotropic lattice); the B/L⁴ term absorbs the leading correction and
    ///     sharpens the estimate of A.
    /// </summary>
    public static (double GInfinity, double A) FitFreeEnergy(IReadOnlyList<int> latticeWidths)
    {
        if (latticeWidths.Count < 3)
        {
            throw new ArgumentException("At least three lattice widths are required for the quadratic fit.", nameof(latticeWidths));
        }

        // Design matrix columns: [1, u, u^2] with u = 1/L^2; solve the 3x3 normal
        // equations directly (small, fixed size — no need for a general solver).
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, t0 = 0, t1 = 0, t2 = 0;
        foreach (var latticeWidth in latticeWidths)
        {
            var u = 1.0 / ((double)latticeWidth * latticeWidth);
            var g = FreeEnergyG(latticeWidth, CriticalCoupling);

            s0 += 1.0;
            s1 += u;
            s2 += u * u;
            s3 += u * u * u;
            s4 += u * u * u * u;
            t0 += g;
            t1 += g * u;
            t2 += g * u * u;
        }

        var (gInfinity, a, _) = SolveSymmetric3x3(
            s0, s1, s2,
            s1, s2, s3,
            s2, s3, s4,
            t0, t1, t2);

        return (gInfinity, a);
    }

    private static (double X, double Y, double Z) SolveSymmetric3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        double b1, double b2, double b3)
    {
        var det = a11 * (a22 * a33 - a23 * a32)
                - a12 * (a21 * a33 - a23 * a31)
                + a13 * (a21 * a32 - a22 * a31);

        var x = (b1 * (a22 * a33 - a23 * a32)
               - a12 * (b2 * a33 - a23 * b3)
               + a13 * (b2 * a32 - a22 * b3)) / det;

        var y = (a11 * (b2 * a33 - a23 * b3)
               - b1 * (a21 * a33 - a23 * a31)
               + a13 * (a21 * b3 - b2 * a31)) / det;

        var z = (a11 * (a22 * b3 - b2 * a32)
               - a12 * (a21 * b3 - b2 * a31)
               + b1 * (a21 * a32 - a22 * a31)) / det;

        return (x, y, z);
    }

    /// <summary>
    ///     Central charge from the finite-size fit: c = 6A/π. For 2D Ising, c → 1/2.
    /// </summary>
    public static double CentralCharge(IReadOnlyList<int> latticeWidths)
    {
        var (_, a) = FitFreeEnergy(latticeWidths);
        return 6.0 * a / Math.PI;
    }

    /// <summary>
    ///     Scaling dimension of the spin field σ at K_c, from the gap between the
    ///     Ramond (periodic) and Neveu–Schwarz (antiperiodic) sector ground states:
    ///     x_σ(L) = (L/4π)·(Σ_NS γ − Σ_R γ). For 2D Ising, x_σ → 1/8.
    /// </summary>
    public static double SigmaDimension(int latticeWidth)
    {
        var sumNs = 0.0;
        var sumR = 0.0;
        for (var r = 0; r < latticeWidth; r++)
        {
            sumNs += Gamma(Math.PI * (2 * r + 1) / latticeWidth, CriticalCoupling);
            sumR += Gamma(2.0 * Math.PI * r / latticeWidth, CriticalCoupling);
        }

        return latticeWidth * (sumNs - sumR) / (4.0 * Math.PI);
    }

    /// <summary>
    ///     Scaling dimension of the energy field ε at K_c, from the lowest
    ///     zero-momentum two-fermion excitation in the Neveu–Schwarz sector (the
    ///     pair at momenta ±π/L): x_ε(L) = (L/π)·γ(π/L). For 2D Ising, x_ε → 1.
    /// </summary>
    public static double EpsilonDimension(int latticeWidth) =>
        latticeWidth * Gamma(Math.PI / latticeWidth, CriticalCoupling) / Math.PI;
}
