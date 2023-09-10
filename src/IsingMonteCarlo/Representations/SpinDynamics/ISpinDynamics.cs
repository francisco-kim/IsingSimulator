namespace IsingMonteCarlo.Representations.SpinDynamics;

public interface ISpinDynamics
{
    void FlipSpin(
        double beta,
        double j,
        double h,
        double? jY);

    void EmptyQueue(
        double beta,
        double j,
        double h,
        double? jY,
        bool verbose = false);
}
