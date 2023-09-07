using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace IsingMonteCarlo.Helpers;

public class DrawHelpers
{
    private const int TwoPowNine = 512;
    private const int ThreePowSix = 729;

    /// <summary>
    /// Creates a bitmap based on data, width, height, stride and pixel format.
    /// </summary>
    /// <param name="sourceData">Byte array of raw source data</param>
    /// <param name="width">Width of the image</param>
    /// <param name="height">Height of the image</param>
    /// <param name="stride">Scanline length inside the data</param>
    /// <param name="pixelFormat">Pixel format</param>
    /// <param name="palette">Color palette</param>
    /// <param name="defaultColor">Default color to fill in on the palette if the given colors don't fully fill it.</param>
    /// <returns>The new image</returns>
    public static Bitmap BuildImage(
        byte[] sourceData,
        int width,
        int height,
        int stride,
        PixelFormat pixelFormat,
        Color[] palette,
        Color? defaultColor)
    {
        var newImage = new Bitmap(width, height, pixelFormat);
        var targetData = newImage.LockBits(
            new Rectangle(x: 0, y: 0, width, height),
            ImageLockMode.WriteOnly,
            newImage.PixelFormat);
        var newDataWidth = ((Image.GetPixelFormatSize(pixelFormat) * width) + 7) / 8;
        // Compensate for possible negative stride on BMP format.
        var isFlipped = stride < 0;
        stride = Math.Abs(stride);
        // Cache these to avoid unnecessary getter calls.
        var targetStride = targetData.Stride;
        var scan0 = targetData.Scan0.ToInt64();
        for (var y = 0; y < height; y++)
        {
            Marshal.Copy(sourceData, y * stride, new IntPtr(scan0 + y * targetStride), newDataWidth);
        }

        newImage.UnlockBits(targetData);
        // Fix negative stride on BMP format.
        if (isFlipped)
        {
            newImage.RotateFlip(RotateFlipType.Rotate180FlipX);
        }

        // For indexed images, set the palette.
        if ((pixelFormat & PixelFormat.Indexed) != 0 && palette is not null)
        {
            var pal = newImage.Palette;
            for (var i = 0; i < pal.Entries.Length; i++)
            {
                if (i < palette.Length)
                {
                    pal.Entries[i] = palette[i];
                }
                else if (defaultColor.HasValue)
                {
                    pal.Entries[i] = defaultColor.Value;
                }
                else
                {
                    break;
                }
            }
            newImage.Palette = pal;
        }
        return newImage;
    }

    public static Bitmap GenerateGrayBitmapFrom2DList(List<int> data)
    {
        // Transform 2-dimensional Int32 array to 1-byte-per-pixel byte array
        var length = data.Count;
        var latticeLength = Convert.ToInt32(Math.Sqrt(length));
        var width = latticeLength;
        var height = width;
        var byteIndex = 0;
        var dataBytes = new byte[length];
        for (var x = 0; x < length; x++)
        {
            // logical AND to be 100% sure the int32 value fits inside
            // the byte even if it contains more data (like, full ARGB).
            dataBytes[byteIndex] = (byte)(((uint)data[x]) & 0xFF);
            // More efficient than multiplying
            byteIndex++;
        }
        // generate palette
        var palette = new Color[256];
        for (var b = 255; b > -1; b--)
        {
            palette[b] = Color.FromArgb(b, b, b);
        }

        // Build image
        return BuildImage(dataBytes, width, height, width, PixelFormat.Format8bppIndexed, palette, null);
    }

    public static Bitmap ResizeBitmap(Bitmap sourceBitmap, int width, int height)
    {
        if (sourceBitmap is null)
        {
            throw new ArgumentNullException(nameof(sourceBitmap));
        }

        if (sourceBitmap.Width < width)
        {
            return ResizeToLargerBitmap(sourceBitmap, width, height);
        }
        else
        {
            return ResizeToSmallerBitmap(sourceBitmap, width, height);
        }
    }

    public static Bitmap ResizeToLargerBitmap(Bitmap sourceBitmap, int width, int height)
    {
        var result = new Bitmap(width, height);
        using var g = Graphics.FromImage(result);
        g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.Half;
        g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
        g.DrawImage(sourceBitmap, 0, 0, width, height);
        return result;
    }

    public static Bitmap ResizeToSmallerBitmap(Bitmap sourceBitmap, int width, int height)
    {
        var result = new Bitmap(width, height);
        using var g = Graphics.FromImage(result);
        // g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
        g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
        g.DrawImage(sourceBitmap, 0, 0, width, height);
        return result;
    }

    public static void SaveBitmapAsPNG(Bitmap bitmap, string filename, bool resize = false, int? size = null)
    {
        var rootDirectory = FileHelpers.GetDataRootDirectory();
        var fullFilename = Path.GetFullPath(Path.Combine(rootDirectory, filename + ".png"));

        if (resize && size is null)
        {
            var resizedBitmap = bitmap.Width % 3 is not 0
                                ? ResizeBitmap(bitmap, TwoPowNine, TwoPowNine)
                                : ResizeBitmap(bitmap, ThreePowSix, ThreePowSix);

            resizedBitmap.Save(fullFilename, ImageFormat.Png);
        }
        else if (resize && size is not null)
        {
            var resizedBitmap = ResizeBitmap(bitmap, (int)size, (int)size);

            resizedBitmap.Save(fullFilename, ImageFormat.Png);
        }
        else
        {
            bitmap.Save(fullFilename, ImageFormat.Png);
        }
    }
}
