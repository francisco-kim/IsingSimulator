namespace IsingMonteCarlo.Representations;

/// <summary>
///     Classical XY model on the 2D square lattice with periodic boundaries:
///     H = -J * sum over nearest-neighbour pairs of cos(theta_i - theta_j),
///     with J = +1 (ferromagnetic). Spins are planar angles theta in [0, 2*pi).
///     Updates are local Metropolis moves with an adaptive proposal window.
///     Site index convention matches the Ising lattice: i = x + y * L.
/// </summary>
public sealed class XYModelSimulator
{
    private const int LatticeSizeLowerBound = 3;
    private const int AcceptanceWindow = 4096;
    private const double TwoPi = 2.0 * Math.PI;

    private readonly double[] _angles;
    private readonly Random _random;

    private double _proposalWidth = Math.PI;
    private int _attemptsInWindow;
    private int _acceptedInWindow;

    public XYModelSimulator(int latticeLength, bool hotStart, int? randomSeed = null)
    {
        if (latticeLength < LatticeSizeLowerBound)
        {
            throw new ArgumentOutOfRangeException(
                nameof(latticeLength),
                $"The lattice length must be at least {LatticeSizeLowerBound}, but {latticeLength} was given.");
        }

        LatticeLength = latticeLength;
        TotalSpinsCount = latticeLength * latticeLength;
        _angles = new double[TotalSpinsCount];
        _random = randomSeed is not null ? new Random((int)randomSeed) : new Random();

        if (hotStart)
        {
            for (var i = 0; i < TotalSpinsCount; i++)
            {
                _angles[i] = TwoPi * _random.NextDouble();
            }
        }
    }

    public int LatticeLength { get; }

    public int TotalSpinsCount { get; }

    public IReadOnlyList<double> Angles => _angles;

    public void RunSteps(long stepCount, double beta)
    {
        var length = LatticeLength;

        for (long step = 0; step < stepCount; step++)
        {
            var site = _random.Next(TotalSpinsCount);
            var x = site % length;
            var y = site - x;

            var right = (x + 1 == length ? 0 : x + 1) + y;
            var left = (x == 0 ? length - 1 : x - 1) + y;
            var down = x + (y + length == TotalSpinsCount ? 0 : y + length);
            var up = x + (y == 0 ? TotalSpinsCount - length : y - length);

            var oldAngle = _angles[site];
            var newAngle = oldAngle + _proposalWidth * (2.0 * _random.NextDouble() - 1.0);

            var deltaEnergy =
                Math.Cos(oldAngle - _angles[right]) - Math.Cos(newAngle - _angles[right])
              + Math.Cos(oldAngle - _angles[left]) - Math.Cos(newAngle - _angles[left])
              + Math.Cos(oldAngle - _angles[down]) - Math.Cos(newAngle - _angles[down])
              + Math.Cos(oldAngle - _angles[up]) - Math.Cos(newAngle - _angles[up]);

            if (deltaEnergy <= 0.0 || _random.NextDouble() <= Math.Exp(-beta * deltaEnergy))
            {
                _angles[site] = newAngle - TwoPi * Math.Floor(newAngle / TwoPi);
                _acceptedInWindow++;
            }

            if (++_attemptsInWindow == AcceptanceWindow)
            {
                AdaptProposalWidth();
            }
        }
    }

    public double GetEnergyPerSite()
    {
        // Each site owns its bonds to the right and downwards, so every bond
        // is counted exactly once.
        var length = LatticeLength;
        var energy = 0.0;

        for (var site = 0; site < TotalSpinsCount; site++)
        {
            var x = site % length;
            var y = site - x;
            var right = (x + 1 == length ? 0 : x + 1) + y;
            var down = x + (y + length == TotalSpinsCount ? 0 : y + length);

            energy -= Math.Cos(_angles[site] - _angles[right])
                    + Math.Cos(_angles[site] - _angles[down]);
        }

        return energy / TotalSpinsCount;
    }

    /// <summary>
    ///     Modulus of the average spin vector, |sum of e^(i*theta)| / N.
    ///     In 2D there is no true long-range order (Mermin-Wagner): on a finite
    ///     lattice this is nonzero below the Kosterlitz-Thouless temperature and
    ///     decays with system size.
    /// </summary>
    public double GetMagnetisationMagnitude()
    {
        var sumCos = 0.0;
        var sumSin = 0.0;
        for (var i = 0; i < TotalSpinsCount; i++)
        {
            sumCos += Math.Cos(_angles[i]);
            sumSin += Math.Sin(_angles[i]);
        }

        return Math.Sqrt(sumCos * sumCos + sumSin * sumSin) / TotalSpinsCount;
    }

    /// <summary>
    ///     Locates vortices (winding +1) and antivortices (winding -1) by
    ///     summing the wrapped angle differences around every elementary
    ///     plaquette. On the torus the two counts balance exactly (the total
    ///     winding vanishes). A plaquette is labelled by its lower-left site.
    /// </summary>
    public (int VortexCount, int AntivortexCount) LocateVortices(
        List<int>? vortexPlaquettes = null,
        List<int>? antivortexPlaquettes = null)
    {
        vortexPlaquettes?.Clear();
        antivortexPlaquettes?.Clear();

        var length = LatticeLength;
        var vortexCount = 0;
        var antivortexCount = 0;

        for (var site = 0; site < TotalSpinsCount; site++)
        {
            var x = site % length;
            var y = site - x;
            var right = (x + 1 == length ? 0 : x + 1) + y;
            var down = x + (y + length == TotalSpinsCount ? 0 : y + length);
            var diagonal = (x + 1 == length ? 0 : x + 1) + (y + length == TotalSpinsCount ? 0 : y + length);

            // Counter-clockwise around the plaquette (site, right, diagonal, down).
            var winding = WrapToPi(_angles[right] - _angles[site])
                        + WrapToPi(_angles[diagonal] - _angles[right])
                        + WrapToPi(_angles[down] - _angles[diagonal])
                        + WrapToPi(_angles[site] - _angles[down]);

            if (winding > Math.PI)
            {
                vortexCount++;
                vortexPlaquettes?.Add(site);
            }
            else if (winding < -Math.PI)
            {
                antivortexCount++;
                antivortexPlaquettes?.Add(site);
            }
        }

        return (vortexCount, antivortexCount);
    }

    private void AdaptProposalWidth()
    {
        var acceptanceRatio = (double)_acceptedInWindow / _attemptsInWindow;
        _proposalWidth = Math.Clamp(
            _proposalWidth * (acceptanceRatio > 0.5 ? 1.1 : 0.9),
            0.05,
            Math.PI);
        _attemptsInWindow = 0;
        _acceptedInWindow = 0;
    }

    private static double WrapToPi(double angleDifference) =>
        angleDifference - TwoPi * Math.Round(angleDifference / TwoPi);
}
