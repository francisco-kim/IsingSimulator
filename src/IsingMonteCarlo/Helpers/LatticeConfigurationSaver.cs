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
        var completePathWithoutFileExtension = Path.GetFullPath(Path.Combine(dataDirectory, filename));

        if (!isBinary)
        {
            File.WriteAllLines(
                completePathWithoutFileExtension + ".dat",
                enumeratedSpins.Select(x => string.Join(separator: ";", x)));
        }
        else
        {
            using var writer = new BinaryWriter(
                new FileStream(path: completePathWithoutFileExtension + ".bin",
                               FileMode.Create));
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
        var filenameData = GetFilenameData(filename);

        temperature = Convert.ToDouble(filenameData[index: 1] + "." + filenameData[index: 2]);
        iterationCount = Convert.ToInt32(filenameData[index: 3]);

        if (filenameData.Last() is "dat")
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

    public static string GetFirstDataFileWithLatticeSizeAndTemperature(int latticeLength, double temperature)
    {
        var dataDirectory = GetDataRootDirectory();
        var firstFileName = new DirectoryInfo(dataDirectory)
                            .EnumerateFiles()
                            .FirstOrDefault(
                                file => filenameContainsCorrectLatticeLengthAndTemperature(
                                    GetFilenameData(file.FullName),
                                    latticeLength,
                                    temperature)) ?? throw new ArgumentNullException(dataDirectory);

        return firstFileName.FullName;
    }

    public static string GetDataRootDirectory()
    {
        var workingDirectory = Directory.GetCurrentDirectory();
        var rootParentDirectory = Directory.GetParent(workingDirectory)
            ?? throw new ArgumentNullException(nameof(workingDirectory));
        var rootDirectory = Directory.GetParent(rootParentDirectory.FullName)
            ?? throw new ArgumentNullException(nameof(rootParentDirectory));

        return Path.GetFullPath(Path.Combine(rootDirectory.FullName, "data"));
    }

    public static (string Filename, string DataDirectory) GetFilename(List<int> spins,
                                                                      double temperature,
                                                                      int iterationCount)
    {
        var latticeLength = Convert.ToInt32(Math.Sqrt(spins.Count));

        return GetFilename(latticeLength, temperature, iterationCount);
    }

    public static (string Filename, string DataDirectory) GetFilename(int latticeLength,
                                                                      double temperature,
                                                                      int iterationCount)
    {
        var dataDirectory = GetDataRootDirectory();

        if (!Directory.Exists(dataDirectory))
        {
            Directory.CreateDirectory(dataDirectory);
        }

        var filename = latticeLength
                     + "_"
                     + $"{temperature:0.0000}"
                     + "_"
                     + iterationCount.ToString(format: "G10", CultureInfo.InvariantCulture);

        return (filename, dataDirectory);
    }

    private static List<string> GetFilenameData(string filename)
    {
        char[] delimiters = { '_', '.' };
        var filenameData = filename.Split(delimiters).ToList();

        if (filenameData.Count != 5)
        {
            throw new ArgumentException(message: "The filename is not formatted correctly.", nameof(filename));
        }

        return filenameData;
    }

    private static bool filenameContainsCorrectLatticeLengthAndTemperature(List<string> filenameData, int latticeSize, double temperature)
    {
        var latticeSizeInFilename = Convert.ToInt32(filenameData[index: 0]);
        var temperatureInFilename = Convert.ToDouble(filenameData[index: 1] + "." + filenameData[index: 2]);

        return latticeSizeInFilename == latticeSize && Math.Abs(temperatureInFilename - temperature) < 1e-5;
    }
}
