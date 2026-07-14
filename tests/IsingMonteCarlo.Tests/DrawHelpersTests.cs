using IsingMonteCarlo.Helpers;

using SixLabors.ImageSharp;

namespace IsingMonteCarlo.Tests;

public class DrawHelpersTests
{
    [Fact]
    public void GenerateGrayBitmapFrom2DList_MapsSpinsToGrayValuesAndEncodesPng()
    {
        var latticeLength = 9;
        var spins = Enumerable.Range(0, latticeLength * latticeLength)
                              .Select(i => i % 2 is 0 ? 1 : -1)
                              .ToList();

        using var image = DrawHelpers.GenerateGrayBitmapFrom2DList(spins);

        Assert.Equal(latticeLength, image.Width);
        Assert.Equal(latticeLength, image.Height);

        // Spin +1 -> gray value 1 (near black); spin -1 -> 255 (white),
        // matching the original System.Drawing palette mapping.
        Assert.Equal(1, image[0, 0].PackedValue);
        Assert.Equal(255, image[1, 0].PackedValue);

        using var stream = new MemoryStream();
        image.SaveAsPng(stream);
        Assert.True(stream.Length > 0);
    }

    [Fact]
    public void ResizeBitmap_ScalesWithNearestNeighbour()
    {
        var spins = Enumerable.Repeat(element: 1, count: 9).ToList();
        using var image = DrawHelpers.GenerateGrayBitmapFrom2DList(spins);

        DrawHelpers.ResizeBitmap(image, width: 81, height: 81);

        Assert.Equal(81, image.Width);
        Assert.Equal(1, image[40, 40].PackedValue);
    }
}
