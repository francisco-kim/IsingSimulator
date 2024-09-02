using System.Linq;
using Avalonia.Controls;
using _2DSimulator.ViewModels;

namespace _2DSimulator.Views;

public partial class MainWindow : Window
{
    public MainWindow()
    {
        InitializeComponent();

        var vm = new MainViewModel();
        DataContext = vm;

        // need to bind this way to use content control on all cells
        Grid.Children.AddRange(
            vm.All.Select(
                spin => new
                    ContentControl { Content = spin }));
    }

    protected override void OnDataContextEndUpdate()
    {
        if (DataContext is MainViewModel vm)
        {
            // only run once the 
            vm.ToggleRun(state: true);
        }
    }
}
