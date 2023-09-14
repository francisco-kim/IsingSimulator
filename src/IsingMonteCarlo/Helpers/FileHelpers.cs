using System.Collections.Concurrent;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.CompilerServices;
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
        var latticeLength = Convert.ToInt32(Math.Sqrt(enumeratedSpins.Count));
        var filename = GetFilename(
            latticeLength,
            temperature,
            iterationCountInMCSweepUnit);
        var dataDirectory = GetDataLatticeLengthSubdirectory(latticeLength);
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
            $"Spins configuration saved in {completePathWithoutFileExtension + extension}.\n");
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
        var dataDirectory = GetDataRootDirectory(new[] { latticeLength.ToString() });

        if (!Directory.Exists(dataDirectory))
        {
            Directory.CreateDirectory(dataDirectory);
        }

        var lastFileName = new DirectoryInfo(dataDirectory)
                            .EnumerateFiles()
                            .Where(file => file.Name.Split('.').Last() == "dat" || file.Name.Split('.').Last() == "bin")
                            .OrderBy(file => file.Name)
                            .LastOrDefault(
                                file => filenameContainsCorrectLatticeLengthAndTemperature(
                                    GetFilenameData(file.Name),
                                    latticeLength,
                                    temperature));

        return lastFileName?.FullName ?? null;
    }

    public static string GetFullPathWithFilename(string filename)
    {
        var filenameData = GetFilenameData(filename);

        var latticeSize = Convert.ToInt32(filenameData[index: 0].Split(Path.DirectorySeparatorChar).Last());

        return GetDataRootDirectory(new[] { Convert.ToString(latticeSize), filename });
    }

    public static string GetDataRootDirectory(IEnumerable<string>? namesToCombine = null)
    {
        var workingDirectory = Directory.GetCurrentDirectory();
        var rootDirectory = Directory.GetParent(workingDirectory).Parent?.Parent?.Parent?.Parent
                         ?? throw new ArgumentNullException(nameof(workingDirectory));

        if (namesToCombine is null)
        {
            return Path.GetFullPath(Path.Combine(rootDirectory.FullName, "data"));
        }

        var enumeratedNamesToCombine = namesToCombine.ToList();
        var combinedPath = new List<string> { rootDirectory.FullName, "data" };
        combinedPath.AddRange(enumeratedNamesToCombine);
        var resultingPath = combinedPath.ToArray() ?? throw new ArgumentNullException(nameof(combinedPath));

        return Path.GetFullPath(Path.Combine(resultingPath));
    }

    public static string GetDataLatticeLengthSubdirectory(int latticeLength)
    {
        // Working directory from launch.json
        var workingDirectory = Directory.GetCurrentDirectory();
        var rootDirectory = Directory.GetParent(workingDirectory)?.Parent?.Parent?.Parent?.Parent
            ?? throw new ArgumentNullException(nameof(workingDirectory));

        return Path.GetFullPath(Path.Combine(rootDirectory.FullName, "data", Convert.ToString(latticeLength)));
    }

    public static string GetFilename(int latticeLength,
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

        return filename;
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

    public static bool FoundFileInLatticeLengthFolder(int latticeLengthToMatchFolderName, out string? filename)
    {
        var dataDirectory = GetDataLatticeLengthSubdirectory(latticeLengthToMatchFolderName);
        var filesInLatticeLengthFolder = Directory
                                         .GetFiles(dataDirectory)
                                         .Where(
                                             file => file.Split(separator: '.').Last() == "dat"
                                                  || file.Split(separator: '.').Last() == "bin");
        foreach (var file in filesInLatticeLengthFolder)
        {
            Console.WriteLine(file);
        }
        Console.WriteLine(value: "Enter file name without the extension (.dat or .bin): ");

        var filenameInput = Console.ReadLine() ?? throw new ArgumentNullException();

        var filenameWithDatExtension = Path.GetFullPath(Path.Combine(dataDirectory, filenameInput + ".dat"));
        var filenameWithBinExtension = Path.GetFullPath(Path.Combine(dataDirectory, filenameInput + ".bin"));
        if (File.Exists(filenameWithBinExtension))
        {
            Console.WriteLine($"Found {filenameWithBinExtension}.\n");
            filename = filenameWithBinExtension;
            return true;
        }

        if (File.Exists(filenameWithDatExtension))
        {
            Console.WriteLine($"Found {filenameWithDatExtension}.\n");
            filename = filenameWithDatExtension;
            return true;
        }

        Console.WriteLine("No such file found.\n");
        filename = null;
        return false;
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
