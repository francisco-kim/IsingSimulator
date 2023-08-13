namespace IsingMonteCarlo.Helpers;

public static class CalculationsHelpers
{
    public static int Modulo(int dividend, int divisor)
    {
        // More efficient than ((dividend % divisor) + divisor) % divisor
        var remainder = dividend % divisor;
        if (remainder < 0)
        {
            remainder += divisor;
        }

        return remainder;
    }
}