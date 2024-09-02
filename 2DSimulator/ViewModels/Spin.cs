using CommunityToolkit.Mvvm.ComponentModel;

namespace _2DSimulator.ViewModels;

public partial class Spin : ViewModelBase
{
    [ObservableProperty, NotifyPropertyChangedFor(nameof(IsSpinUp))]
    private int spinValue;

    public Spin(int x, int y, int spinValue)
    {
        X = x;
        Y = y;
        SpinValue = spinValue;
    }

    public int X { get; }

    public int Y { get; }

    public bool IsSpinUp => SpinValue == +1;
}
