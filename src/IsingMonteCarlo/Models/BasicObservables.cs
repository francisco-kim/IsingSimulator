namespace IsingMonteCarlo.Models;

public record class BasicObservables(
    double Magnetisation,
    double MagnetisationSquared,
    double MagnetisationAbsolute,
    double Energy,
    double CorrelationLengthX,
    double CorrelationLengthY,
    double RenormalisedCorrelationLength,
    double Susceptibility,
    double SpecificHeat,
    List<double> MagnetisationList,
    List<double> MagnetisationSquaredList,
    List<double> MagnetisationAbsoluteList,
    List<double> EnergyList);
