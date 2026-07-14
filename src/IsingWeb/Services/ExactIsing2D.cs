namespace IsingWeb.Services;

/// <summary>
///     Exact results for the square-lattice Ising ferromagnet with |J| = 1, h = 0.
/// </summary>
public static class ExactIsing2D
{
    /// <summary>T_c = 2 / ln(1 + sqrt(2)) ≈ 2.2692.</summary>
    public static readonly double CriticalTemperature = 2.0 / Math.Log(1.0 + Math.Sqrt(2.0));

    /// <summary>
    ///     Onsager–Yang spontaneous magnetisation:
    ///     M(T) = (1 - sinh(2/T)^-4)^(1/8) for T &lt; T_c, 0 above.
    /// </summary>
    public static double Magnetisation(double temperature)
    {
        if (temperature <= 0.0)
        {
            return 1.0;
        }

        if (temperature >= CriticalTemperature)
        {
            return 0.0;
        }

        var sinh = Math.Sinh(2.0 / temperature);
        return Math.Pow(1.0 - Math.Pow(sinh, -4.0), 0.125);
    }
}
