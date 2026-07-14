namespace IsingWeb.Services;

public enum LatticeColourMode
{
    Spins,
    Domains
}

/// <summary>
///     Converts a spin configuration into an RGBA pixel buffer.
///     Rendering convention: spin index i sits at pixel (i % L, i / L).
/// </summary>
public sealed class LatticeRenderer
{
    // Ink scheme: spin up near-black, spin down near-white — matches the
    // repository's classic grayscale PNGs and reads on light and dark pages.
    public static readonly (byte R, byte G, byte B) SpinUpColour = (0x16, 0x18, 0x1d);
    public static readonly (byte R, byte G, byte B) SpinDownColour = (0xfc, 0xfc, 0xfb);

    // Orange (categorical slot 8) — distinct from both spin colours.
    public static readonly (byte R, byte G, byte B) ClusterColour = (0xeb, 0x68, 0x34);

    private byte[] _rgba = Array.Empty<byte>();
    private int[] _domainRoot = Array.Empty<int>();

    public byte[] Buffer => _rgba;

    public void EnsurePixelCount(int pixelCount)
    {
        if (_rgba.Length != pixelCount * 4)
        {
            _rgba = new byte[pixelCount * 4];
        }
    }

    public void RenderSpins(IReadOnlyList<int> spins)
    {
        EnsurePixelCount(spins.Count);
        for (var i = 0; i < spins.Count; i++)
        {
            var (r, g, b) = spins[i] > 0 ? SpinUpColour : SpinDownColour;
            WritePixel(i, r, g, b);
        }
    }

    /// <summary>
    ///     Renders planar spins coloured by their phase: hue = angle. Used by
    ///     the XY model page.
    /// </summary>
    public void RenderPhases(IReadOnlyList<double> angles)
    {
        EnsurePixelCount(angles.Count);
        for (var i = 0; i < angles.Count; i++)
        {
            var hue = angles[i] * (180.0 / Math.PI);
            var (r, g, b) = HslToRgb(hue, saturation: 0.75, lightness: 0.55);
            WritePixel(i, r, g, b);
        }
    }

    /// <summary>
    ///     Marks plaquettes (labelled by their lower-left site) with a 5-pixel
    ///     cross — white for vortices, near-black for antivortices. Call after
    ///     <see cref="RenderPhases" />.
    /// </summary>
    public void MarkPlaquettes(
        IReadOnlyList<int> plaquetteSites,
        int latticeLength,
        bool white,
        bool drawCross = true)
    {
        var height = _rgba.Length / 4 / latticeLength;
        var (r, g, b) = white ? ((byte)255, (byte)255, (byte)255) : ((byte)10, (byte)10, (byte)10);

        foreach (var site in plaquetteSites)
        {
            WritePixel(site, r, g, b);

            if (!drawCross)
            {
                continue;
            }

            var x = site % latticeLength;
            var y = site / latticeLength;
            WritePixel((x + 1) % latticeLength + y * latticeLength, r, g, b);
            WritePixel((x + latticeLength - 1) % latticeLength + y * latticeLength, r, g, b);
            WritePixel(x + (y + 1) % height * latticeLength, r, g, b);
            WritePixel(x + (y + height - 1) % height * latticeLength, r, g, b);
        }
    }

    /// <summary>
    ///     Overlays the given sites (e.g. the currently growing Wolff cluster)
    ///     in the cluster accent colour. Call after <see cref="RenderSpins" />.
    /// </summary>
    public void HighlightSites(IReadOnlyList<int> sites)
    {
        var pixelCount = _rgba.Length / 4;
        foreach (var site in sites)
        {
            if (site >= 0 && site < pixelCount)
            {
                WritePixel(site, ClusterColour.R, ClusterColour.G, ClusterColour.B);
            }
        }
    }

    /// <summary>
    ///     Colours each connected same-spin domain (periodic boundaries) with its
    ///     own hue; up-domains are darker, down-domains lighter, so the spin sign
    ///     stays readable.
    /// </summary>
    public void RenderDomains(IReadOnlyList<int> spins, int latticeLength)
    {
        var count = spins.Count;
        EnsurePixelCount(count);

        if (_domainRoot.Length != count)
        {
            _domainRoot = new int[count];
        }

        for (var i = 0; i < count; i++)
        {
            _domainRoot[i] = i;
        }

        var height = count / latticeLength;
        for (var i = 0; i < count; i++)
        {
            var x = i % latticeLength;
            var y = i / latticeLength;
            var right = (x + 1) % latticeLength + y * latticeLength;
            var down = x + ((y + 1) % height) * latticeLength;

            if (spins[i] == spins[right])
            {
                Union(i, right);
            }

            if (height > 1 && spins[i] == spins[down])
            {
                Union(i, down);
            }
        }

        for (var i = 0; i < count; i++)
        {
            var root = Find(i);

            // Golden-angle hue walk keyed on the domain root: nearby roots get
            // well-separated hues. Lightness encodes the spin sign.
            var hue = root * 137.508 % 360.0;
            var lightness = spins[i] > 0 ? 0.38 : 0.72;
            var (r, g, b) = HslToRgb(hue, saturation: 0.55, lightness);
            WritePixel(i, r, g, b);
        }
    }

    private void WritePixel(int index, byte r, byte g, byte b)
    {
        var offset = index * 4;
        _rgba[offset] = r;
        _rgba[offset + 1] = g;
        _rgba[offset + 2] = b;
        _rgba[offset + 3] = 255;
    }

    private int Find(int site)
    {
        var root = site;
        while (_domainRoot[root] != root)
        {
            root = _domainRoot[root];
        }

        while (_domainRoot[site] != root)
        {
            (site, _domainRoot[site]) = (_domainRoot[site], root);
        }

        return root;
    }

    private void Union(int a, int b)
    {
        var rootA = Find(a);
        var rootB = Find(b);
        if (rootA != rootB)
        {
            _domainRoot[rootB] = rootA;
        }
    }

    private static (byte R, byte G, byte B) HslToRgb(double hue, double saturation, double lightness)
    {
        var c = (1.0 - Math.Abs(2.0 * lightness - 1.0)) * saturation;
        var hPrime = hue / 60.0;
        var x = c * (1.0 - Math.Abs(hPrime % 2.0 - 1.0));
        var (r1, g1, b1) = (int)hPrime switch
        {
            0 => (c, x, 0.0),
            1 => (x, c, 0.0),
            2 => (0.0, c, x),
            3 => (0.0, x, c),
            4 => (x, 0.0, c),
            _ => (c, 0.0, x)
        };
        var m = lightness - c / 2.0;

        return ((byte)((r1 + m) * 255.0), (byte)((g1 + m) * 255.0), (byte)((b1 + m) * 255.0));
    }
}
