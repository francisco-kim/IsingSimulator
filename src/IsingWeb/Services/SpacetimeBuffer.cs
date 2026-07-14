namespace IsingWeb.Services;

/// <summary>
///     Scrolling spacetime diagram for the 1D chain: each appended row is the
///     chain at one Monte Carlo sweep; rows are kept in a circular buffer and
///     copied out oldest-first (top = past, bottom = present).
/// </summary>
public sealed class SpacetimeBuffer
{
    private byte[] _rows = Array.Empty<byte>();
    private byte[] _ordered = Array.Empty<byte>();
    private int _width;
    private int _head;
    private int _filledRows;

    public int Width => _width;

    public int Height { get; private set; }

    public void Configure(int width, int height)
    {
        _width = width;
        Height = height;
        _rows = new byte[width * height * 4];
        _ordered = new byte[width * height * 4];
        _head = 0;
        _filledRows = 0;

        // Start on the surface colour rather than black.
        var (r, g, b) = LatticeRenderer.SpinDownColour;
        for (var i = 0; i < width * height; i++)
        {
            var offset = i * 4;
            _rows[offset] = r;
            _rows[offset + 1] = g;
            _rows[offset + 2] = b;
            _rows[offset + 3] = 255;
        }
    }

    public void AppendRow(IReadOnlyList<int> spins)
    {
        if (spins.Count != _width)
        {
            return;
        }

        var offset = _head * _width * 4;
        for (var x = 0; x < _width; x++)
        {
            var (r, g, b) = spins[x] > 0 ? LatticeRenderer.SpinUpColour : LatticeRenderer.SpinDownColour;
            _rows[offset] = r;
            _rows[offset + 1] = g;
            _rows[offset + 2] = b;
            _rows[offset + 3] = 255;
            offset += 4;
        }

        _head = (_head + 1) % Height;
        _filledRows = Math.Min(_filledRows + 1, Height);
    }

    /// <summary>Returns the buffer with rows ordered oldest-first, ready to blit.</summary>
    public byte[] OrderedBuffer()
    {
        var rowBytes = _width * 4;
        var oldest = _filledRows < Height ? 0 : _head;

        for (var row = 0; row < Height; row++)
        {
            var source = (oldest + row) % Height;
            Array.Copy(_rows, source * rowBytes, _ordered, row * rowBytes, rowBytes);
        }

        return _ordered;
    }
}
