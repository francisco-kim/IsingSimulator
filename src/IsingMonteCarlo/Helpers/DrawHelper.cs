using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace IsingMonteCarlo.Helpers
{
    public class DrawHelper
    {
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
        public static Bitmap BuildImage(Byte[] sourceData, Int32 width, Int32 height, Int32 stride, PixelFormat pixelFormat, Color[] palette, Color? defaultColor)
        {
            Bitmap newImage = new Bitmap(width, height, pixelFormat);
            BitmapData targetData = newImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, newImage.PixelFormat);
            Int32 newDataWidth = ((Image.GetPixelFormatSize(pixelFormat) * width) + 7) / 8;
            // Compensate for possible negative stride on BMP format.
            Boolean isFlipped = stride < 0;
            stride = Math.Abs(stride);
            // Cache these to avoid unnecessary getter calls.
            Int32 targetStride = targetData.Stride;
            Int64 scan0 = targetData.Scan0.ToInt64();
            for (Int32 y = 0; y < height; y++)
                Marshal.Copy(sourceData, y * stride, new IntPtr(scan0 + y * targetStride), newDataWidth);
            newImage.UnlockBits(targetData);
            // Fix negative stride on BMP format.
            if (isFlipped)
                newImage.RotateFlip(RotateFlipType.Rotate180FlipX);
            // For indexed images, set the palette.
            if ((pixelFormat & PixelFormat.Indexed) != 0 && palette != null)
            {
                ColorPalette pal = newImage.Palette;
                for (Int32 i = 0; i < pal.Entries.Length; i++)
                {
                    if (i < palette.Length)
                        pal.Entries[i] = palette[i];
                    else if (defaultColor.HasValue)
                        pal.Entries[i] = defaultColor.Value;
                    else
                        break;
                }
                newImage.Palette = pal;
            }
            return newImage;
        }

        public static Bitmap FromTwoDimIntArrayGray(List<int> data)
        {
            // Transform 2-dimensional Int32 array to 1-byte-per-pixel byte array
            var length = data.Count;
            var latticeLength = Convert.ToInt32(Math.Sqrt(length));
            var width = latticeLength;
            var height = width;
            int byteIndex = 0;
            byte[] dataBytes = new byte[length];
            for (int x = 0; x < length; x++)
            {
                // logical AND to be 100% sure the int32 value fits inside
                // the byte even if it contains more data (like, full ARGB).
                dataBytes[byteIndex] = (byte)(((uint)data[x]) & 0xFF);
                // More efficient than multiplying
                byteIndex++;
            }
            // generate palette
            var palette = new Color[256];
            for (int b = 255; b > -1; b--)
                palette[b] = Color.FromArgb(b, b, b);
            // Build image
            return BuildImage(dataBytes, width, height, width, PixelFormat.Format8bppIndexed, palette, null);
        }

        public static Bitmap ResizeToLargerBitmap(Bitmap sourceBitmap, int width, int height)
        {
            Bitmap result = new Bitmap(width, height);
            using (Graphics g = Graphics.FromImage(result))
            {
                g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.Half;
                g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
                g.DrawImage(sourceBitmap, 0, 0, width, height);
            }
            return result;
        }

        public static Bitmap ResizeToSmallerBitmap(Bitmap sourceBitmap, int width, int height)
        {
            Bitmap result = new Bitmap(width, height);
            using (Graphics g = Graphics.FromImage(result))
            {
                g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
                g.DrawImage(sourceBitmap, 0, 0, width, height);
            }
            return result;
        }

        public static void SaveBmpAsPNG(Bitmap bitmap, string filename)
        {
            var rootDirectory = LatticeConfigurationSaver.GetDataRootDirectory();
            var fullFilename = Path.GetFullPath(Path.Combine(rootDirectory, filename + ".png"));

            bitmap.Save(fullFilename, ImageFormat.Png);
        }
    }
}