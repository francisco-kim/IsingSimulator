using SixLabors.ImageSharp;
using SixLabors.ImageSharp.PixelFormats;
using SixLabors.ImageSharp.Processing;

namespace IsingMonteCarlo.Helpers;

public static class DrawHelpers
{
    private const int TwoPowNine = 512;
    private const int ThreePowSix = 729;

    /// <summary>
    ///     Renders a square spin configuration as an 8-bit grayscale image,
    ///     one pixel per spin: +1 maps to gray value 1 (near black) and -1 to 255 (white).
    /// </summary>
    public static Image<L8> GenerateGrayBitmapFrom2DList(List<int> data)
    {
        var latticeLength = Convert.ToInt32(Math.Sqrt(data.Count));

        if (latticeLength * latticeLength != data.Count)
        {
            throw new ArgumentException(
                $"The spin configuration must be square, but {data.Count} spins were given.",
                nameof(data));
        }

        var image = new Image<L8>(latticeLength, latticeLength);
        image.ProcessPixelRows(accessor =>
        {
            for (var y = 0; y < latticeLength; y++)
            {
                var row = accessor.GetRowSpan(y);
                for (var x = 0; x < latticeLength; x++)
                {
                    row[x] = new L8((byte)(((uint)data[y * latticeLength + x]) & 0xFF));
                }
            }
        });

        return image;
    }

    public static Image<L8> ResizeBitmap(Image<L8> sourceImage, int width = TwoPowNine, int height = TwoPowNine)
    {
        if (sourceImage is null)
        {
            throw new ArgumentNullException(nameof(sourceImage));
        }

        // Nearest-neighbour keeps the spins as crisp blocks in both directions.
        sourceImage.Mutate(context => context.Resize(new ResizeOptions
        {
            Size = new Size(width, height),
            Sampler = KnownResamplers.NearestNeighbor,
            Mode = ResizeMode.Stretch
        }));

        return sourceImage;
    }

    public static void SaveBitmapAsPNGInSpecifiedFolder(
        Image<L8> image,
        int latticeLength,
        double temperature,
        IEnumerable<string>? subPathNames = null,
        bool resize = false,
        int? size = null)
    {
        var filename = FileHelpers.GetFilename(latticeLength, temperature, iterationCountInMCSweepUnit: 0);

        var imagesDirectory =
            FileHelpers.GetDataRootDirectory(subPathNames ?? new[] { Convert.ToString(latticeLength), "images" });

        if (!Directory.Exists(imagesDirectory))
        {
            Directory.CreateDirectory(imagesDirectory);
        }

        var fullFilename = Path.GetFullPath(Path.Combine(imagesDirectory, filename + ".png"));

        if (File.Exists(fullFilename))
        {
            var iterationCount = 0;
            while (File.Exists(fullFilename))
            {
                ++iterationCount;
                filename = FileHelpers.GetFilename(latticeLength, temperature, iterationCount);
                fullFilename = Path.GetFullPath(Path.Combine(imagesDirectory, filename + ".png"));
            }
        }

        if (resize && size is null)
        {
            var targetSize = image.Width % 3 is not 0 ? TwoPowNine : ThreePowSix;
            ResizeBitmap(image, targetSize, targetSize);
        }
        else if (resize && size is not null)
        {
            ResizeBitmap(image, (int)size, (int)size);
        }

        image.SaveAsPng(fullFilename);

        Console.WriteLine($"File {fullFilename} saved.");
    }

    public static void SaveBitmapAsPNG(
        Image<L8> image,
        int latticeLength,
        double temperature,
        bool resize = false,
        int? size = null) =>
        SaveBitmapAsPNGInSpecifiedFolder(
            image,
            latticeLength,
            temperature,
            resize: resize,
            size: size);
}
