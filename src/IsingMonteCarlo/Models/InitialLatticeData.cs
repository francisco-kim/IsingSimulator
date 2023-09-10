namespace IsingMonteCarlo.Models;

public record class InitialLatticeData(List<int> SpinConfiguration, double boltzmannTemperature, int previousIterationCount)
{
    public static implicit operator (List<int> SpinConfiguration, double boltzmannTemperature, int previousIterationCount)(InitialLatticeData value)
    {
        return (value.SpinConfiguration, value.boltzmannTemperature, value.previousIterationCount);
    }

    public static implicit operator InitialLatticeData((List<int> SpinConfiguration, double boltzmannTemperature, int previousIterationCount) value)
    {
        return new InitialLatticeData(value.SpinConfiguration, value.boltzmannTemperature, value.previousIterationCount);
    }
}