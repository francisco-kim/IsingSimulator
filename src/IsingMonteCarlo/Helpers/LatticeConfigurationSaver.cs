using System.Globalization;

namespace IsingMonteCarlo.Helpers;

public static class LatticeConfigurationSaver
{
    public static void SaveLattice(
        IEnumerable<int> spins,
        double temperature,
        int iterationCount,
        bool isBinary = false)
    {
        var enumeratedSpins = spins.ToList() ?? throw new ArgumentException(nameof(spins));
        var (filename, dataDirectory) = GetFilename(enumeratedSpins, temperature, iterationCount);

        if (!isBinary)
        {
            File.WriteAllLines(
                dataDirectory + "\\" + filename + ".dat",
                enumeratedSpins.Select(x => string.Join(separator: ";", x)));
        }
        else
        {
            using var writer = new BinaryWriter(
                new FileStream(path: dataDirectory + "\\" + filename + ".bin", FileMode.Create));
            foreach (var item in enumeratedSpins)
            {
                writer.Write(item);
            }
        }
    }

    public static List<int> LoadLattice(
        string filename,
        out double temperature,
        out int iterationCount)
    {
        char[] delimiters = { '_', '.' };
        var filenameData = filename.Split(delimiters).ToList();

        if (filenameData.Count != 5)
        {
            throw new ArgumentException(message: "The filename is not formatted correctly.", nameof(filename));
        }

        temperature = Convert.ToDouble(filenameData[index: 1] + "." +filenameData[index: 2]);
        iterationCount = Convert.ToInt32(filenameData[index: 3]);

        if (filenameData.Last() is ".dat")
        {
            return File.ReadLines(filename)
                       .SelectMany(x => x.Split(separator: ',').Select(int.Parse))
                       .ToList();
        }

        var result = new List<int>();

        using var reader = new BinaryReader(new FileStream(filename, FileMode.Open));
        while (reader.BaseStream.Position < reader.BaseStream.Length)
        {
            result.Add(reader.ReadInt32());
        }

        return result;
    }

    private static (string Filename, string DataDirectory) GetFilename(List<int> spins, double temperature, int iterationCount)
    {
        var workingDirectory = Environment.CurrentDirectory;
        var dataDirectory = Directory.GetParent(workingDirectory).Parent.Parent.Parent.Parent.FullName + "\\data";
        Directory.CreateDirectory(dataDirectory);

        var filename = Math.Sqrt(spins.Count).ToString(format: "G10", CultureInfo.InvariantCulture)
                     + "_"
                     + $"{temperature:0.0000}"
                     + "_"
                     + iterationCount;

        return (filename, dataDirectory);
    }
}
