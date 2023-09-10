using IsingMonteCarlo.Representations;

namespace IsingMonteCarlo.Helpers;

public static class FileHelpers
{
    public static void SaveSpinConfiguration(
        IEnumerable<int> spins,
        double temperature,
        long iterationCountInMCSweepUnit,
        bool isInByte = true)
    {
        var enumeratedSpins = spins.ToList() ?? throw new ArgumentException(nameof(spins));
        var (filename, dataDirectory) = GetFilename(
            Convert.ToInt32(Math.Sqrt(enumeratedSpins.Count)),
            temperature,
            iterationCountInMCSweepUnit);
        var completePathWithoutFileExtension = Path.GetFullPath(Path.Combine(dataDirectory, filename));

        var extension = ".dat";
        if (!isInByte)
        {
            File.WriteAllLines(
                completePathWithoutFileExtension + extension,
                enumeratedSpins.Select(x => string.Join(separator: ";", x)));
        }
        else
        {
            extension = ".bin";
            var spinsInByte = enumeratedSpins.Select(s => s == 1 ? Convert.ToByte(true) : Convert.ToByte(false))
                                             .ToArray();
            using var fs = new FileStream(completePathWithoutFileExtension + extension, FileMode.Create, FileAccess.Write);
            fs.Write(spinsInByte, 0, spinsInByte.Length);

            //using var writer = new BinaryWriter(
            //    new FileStream(path: completePathWithoutFileExtension + extension,
            //                   FileMode.Create));
            //foreach (var item in enumeratedSpins)
            //{
            //    writer.Write(item);
            //}
        }

        Console.WriteLine(
            $"Filename: {completePathWithoutFileExtension + extension}\n");
    }

    public static List<int> LoadSpinConfiguration(
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

        var fileStream = new FileStream(filename, FileMode.Open, FileAccess.Read);
        var reader = new BinaryReader(fileStream);
        var bytesCount = new FileInfo(filename).Length;
        var result = reader.ReadBytes((int)bytesCount);
        return result.Select(s => s == 1 ? 1 : -1).ToList();

        //var result = new List<int>();

        //using var reader = new BinaryReader(new FileStream(filename, FileMode.Open));
        //while (reader.BaseStream.Position < reader.BaseStream.Length)
        //{
        //    result.Add(reader.ReadBytes() == true ? 1 : -1);
        //}

        //return result;
    }

    public static NearestNeighbourNDIsingLattice<int> LoadLattice(string filename, int dimension)
    {
        var spinConfiguration = LoadSpinConfiguration(GetFullPathWithFilename(filename), out _, out _);
        var latticeLength = Convert.ToInt32(Math.Pow(spinConfiguration.Count, 1.0 / dimension));

        return new NearestNeighbourNDIsingLattice<int>(dimension, latticeLength, spinConfiguration);
    }

    public static string? GetLastDataFileWithLatticeSizeAndTemperature(int latticeLength, double temperature)
    {
        var dataDirectory = GetDataRootDirectory();
        var firstFileName = new DirectoryInfo(dataDirectory)
                            .EnumerateFiles()
                            .Where(file => file.Name.Split('.').Last() == "dat" || file.Name.Split('.').Last() == "bin")
                            .LastOrDefault(
                                file => filenameContainsCorrectLatticeLengthAndTemperature(
                                    GetFilenameData(file.Name),
                                    latticeLength,
                                    temperature));

        return firstFileName?.FullName;
    }

    public static string GetFullPathWithFilename(string filename)
    {
        var filenameData = GetFilenameData(filename);

        var latticeSize = Convert.ToInt32(filenameData[index: 0]);

        return GetDataRootDirectory(new string[] { Convert.ToString(latticeSize), filename });
    }

    public static string GetDataRootDirectory(string[]? namesToCombine = null)
    {
        var workingDirectory = Directory.GetCurrentDirectory();
        var rootParentDirectory = Directory.GetParent(workingDirectory)
            ?? throw new ArgumentNullException(nameof(workingDirectory));
        var rootDirectory = Directory.GetParent(rootParentDirectory.FullName)
            ?? throw new ArgumentNullException(nameof(rootParentDirectory));

        if (namesToCombine is null)
        {
            return Path.GetFullPath(Path.Combine(rootDirectory.FullName, "data"));
        }
        
        var combinedPath = new[] {rootDirectory.FullName, "data"};
        foreach (var name in namesToCombine)
        {
            combinedPath.Append(name);
        }

        return Path.GetFullPath(Path.Combine(combinedPath));
    }

    public static string GetDataLatticeLengthSubdirectory(int latticeLength)
    {
        var workingDirectory = Directory.GetCurrentDirectory();
        var rootParentDirectory = Directory.GetParent(workingDirectory)
            ?? throw new ArgumentNullException(nameof(workingDirectory));
        var rootDirectory = Directory.GetParent(rootParentDirectory.FullName)
            ?? throw new ArgumentNullException(nameof(rootParentDirectory));

        return Path.GetFullPath(Path.Combine(rootDirectory.FullName, "data", Convert.ToString(latticeLength)));
    }

    public static (string Filename, string DataDirectory) GetFilename(int latticeLength,
                                                                      double temperature,
                                                                      long iterationCountInMCSweepUnit)
    {
        var dataDirectory = GetDataLatticeLengthSubdirectory(latticeLength);

        if (!Directory.Exists(dataDirectory))
        {
            Directory.CreateDirectory(dataDirectory);
        }

        string filename;
        if (iterationCountInMCSweepUnit < 0)
        {
            filename = latticeLength
                         + "_"
                         + $"{temperature:0.00000}";
        }
        else
        {
            filename = latticeLength
                         + "_"
                         + $"{temperature:0.00000}"
                         + "_"
                         //  + iterationCount.ToString(format: "G10", CultureInfo.InvariantCulture);
                         + iterationCountInMCSweepUnit.ToString(format: "D8");
        }

        return (filename, dataDirectory);
    }

    public static List<string> GetFilenameData(string filename)
    {
        char[] delimiters = { '_', '.' };
        var filenameData = filename.Split(delimiters).ToList();

        if (filenameData.Count != 5)
        {
            throw new ArgumentException(message: "The filename is not formatted correctly.", nameof(filename));
        }

        return filenameData;
    }

    private static bool filenameContainsCorrectLatticeLengthAndTemperature(
        List<string> filenameData,
        int latticeSize,
        double temperature)
    {
        var latticeSizeInFilename = Convert.ToInt32(filenameData[index: 0]);
        var temperatureInFilename = Convert.ToDouble(filenameData[index: 1] + "." + filenameData[index: 2]);

        return latticeSizeInFilename == latticeSize && Math.Abs(temperatureInFilename - temperature) < 1e-5;
    }
}
