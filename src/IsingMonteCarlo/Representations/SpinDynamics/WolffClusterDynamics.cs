namespace IsingMonteCarlo.Representations.SpinDynamics;

public sealed class WolffClusterDynamics : ISpinDynamics
{
    private readonly IHamiltonian<int> _hamiltonian;
    private readonly int _totalSpinsCount;
    private readonly Queue<int> _clusterQueue;
    private readonly List<int> _currentClusterSites;
    private readonly Random _random;

    private int _referenceSpin;

    public WolffClusterDynamics(IHamiltonian<int> hamiltonian, int? randomSeed = null)
    {
        _hamiltonian = hamiltonian ?? throw new ArgumentNullException(nameof(hamiltonian));
        _totalSpinsCount = hamiltonian.Lattice.TotalSpinsCount;
        _clusterQueue = new Queue<int>();
        _currentClusterSites = new List<int>();
        _random = randomSeed is not null ? new Random((int)randomSeed) : new Random();
    }

    /// <summary>
    ///     True while a cluster is still being grown, i.e. sites remain in the queue.
    /// </summary>
    public bool HasPendingCluster => _clusterQueue.Count > 0;

    /// <summary>
    ///     The sites flipped as part of the cluster currently being grown
    ///     (cleared when a new cluster is seeded). Intended for visualisation.
    /// </summary>
    public IReadOnlyList<int> CurrentClusterSites => _currentClusterSites;

    public void FlipSpin(
        double beta,
        double j,
        double h,
        double? jY)
    {
        if (h is not 0.0)
        {
            throw new ArgumentException(
                "The Wolff single-cluster algorithm is not allowed "
              + $"if the external field is present ({nameof(h)} = {h} here).",
                nameof(h));
        }

        // If a cluster has not been chosen yet
        if (_clusterQueue.Count is 0)
        {
            var chosenSite = _random.Next(minValue: 0, _totalSpinsCount);
            _referenceSpin = _hamiltonian.Lattice.Spins[chosenSite];
            _clusterQueue.Enqueue(chosenSite);
            _currentClusterSites.Clear();
            _currentClusterSites.Add(chosenSite);

            // Flip reference spin of the cluster
            _hamiltonian.FlipSpinWithPropertiesUpdate(chosenSite, j, h, jY);
        }

        var siteIndexToConsider = _clusterQueue.Dequeue();
        var neighboursIndices = _hamiltonian.Lattice.NeighboursIndices[siteIndexToConsider];

        if (jY is null)
        {
            foreach (var neighbour in neighboursIndices)
            {
                if (_hamiltonian.Lattice.Spins[neighbour] != _referenceSpin)
                {
                    continue;
                }

                var addToCluster = _random.NextDouble() <= 1.0 - Math.Exp(2.0 * beta * j);
                if (addToCluster)
                {
                    _clusterQueue.Enqueue(neighbour);
                    _currentClusterSites.Add(neighbour);
                    _hamiltonian.FlipSpinWithPropertiesUpdate(neighbour, j, h, jY);
                }
            }
        }
        else
        {
            var jYCoupling = (double)jY;
            var couplingsList = new List<double> { j, j, jYCoupling, jYCoupling };
            var neighboursWithCouplings = neighboursIndices.Zip(couplingsList);

            foreach (var (neighbour, couplingStrength) in neighboursWithCouplings)
            {
                if (_hamiltonian.Lattice.Spins[neighbour] != _referenceSpin)
                {
                    continue;
                }

                var addToCluster = _random.NextDouble() <= 1.0 - Math.Exp(2.0 * beta * couplingStrength);
                if (addToCluster)
                {
                    _clusterQueue.Enqueue(neighbour);
                    _currentClusterSites.Add(neighbour);
                    _hamiltonian.FlipSpinWithPropertiesUpdate(neighbour, j, h, jY);
                }
            }

            // var flipNeighbouringSpins = neighboursIndices.Take(2)
            //     .Select(s => random.NextDouble() <= 1 - Math.Exp(-2.0 * _beta * j))
            //     .ToList();
            // flipNeighbouringSpins.AddRange(neighboursIndices.Skip(2)
            //     .Select(s => random.NextDouble() <= 1 - Math.Exp(-2.0 * _beta * (double)jY))
            //     .ToList());
        }
    }

    public void EmptyQueue(
        double beta,
        double j,
        double h,
        double? jY,
        bool verbose = false)
    {
        if (_clusterQueue.Count != 0)
        {
            if (verbose)
            {
                Console.Write($"Dequeuing last cluster: {_clusterQueue.Count} spins remaining...");
                Console.Write("\r");
            }
            while (_clusterQueue.Count > 0)
            {
                FlipSpin(beta, j, h, jY);
            }
        }
    }
}
