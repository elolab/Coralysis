# Coralysis 0.99.10

- Main changes: 

    - `ReferenceMapping()`: implementing the `label.prune.cutoff` parameter to prune low-confidence predicted cell labels based on the confidence probability scores
    

# Coralysis 0.99.4

- Main changes: 

    - `RunParallelDivisiveICP()`: adopting `BiocParallel` interface for parallelization
    
    - `RunParallelDivisiveICP()`: added the `RNGseed` parameter (set to `123` by default) to ensure reproducibility during parallelization via `BiocParallel`
    

# Coralysis 0.99.0

First release 🎉

- Functions:
    - `AggregateDataByBatch()`
    - `RunParallelDivisiveICP()`
    - `ReferenceMapping()`
    - `GetCellClusterProbability()`
    - `SummariseCellClusterProbability()`
    - `CellClusterProbabilityDistribution()`
    - `BinCellClusterProbability()`
    - `TabulateCellBinsByGroup()`
    - `CellBinsFeatureCorrelation()`
    - `FindClusterMarkers()`
    - `FindAllClusterMarkers()`
    - `MajorityVotingFeatures()`
    - `RunPCA()`
    - `RunTSNE()`
    - `RunUMAP()`
    - `PCAElbowPlot()`
    - `PlotDimRed()`
    - `PlotExpression()`
    - `PlotClusterTree()`
    - `VlnPlot()`
    - `HeatmapFeatures()`
    - `PrepareData()`
    - `GetFeatureCoefficients()`
- Vignettes:
    - [Integration](https://elolab.github.io/Coralysis/articles/01_Integration.html)
    - [Reference-mapping](https://elolab.github.io/Coralysis/articles/02_RefMap.html)
    - [Cell States](https://elolab.github.io/Coralysis/articles/CellState.html)
