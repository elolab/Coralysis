#' @export
setGeneric("PrepareData",signature = "object",
           function(object) {
             standardGeneric("PrepareData")
           })

#' @export
setGeneric("AggregateDataByBatch", signature = "object",
           function(object, batch.label, 
                    nhvg = 2000L, p = 30L, 
                    ...) {
               standardGeneric("AggregateDataByBatch")
           })

#' @export
setGeneric("RunPCA", signature = "object",
           function(object,
                    p = 50,
                    scale = TRUE,
                    center = TRUE,
                    threshold = 0,
                    method = "irlba", 
                    return.model = FALSE, 
                    select.icp.tables = NULL) {
             standardGeneric("RunPCA")
           })

#' @export
setGeneric("PCAElbowPlot", signature = "object",
           function(object,
                    return.plot=FALSE) {
             standardGeneric("PCAElbowPlot")
           })

#' @export
setGeneric("RunUMAP", signature = "object",
           function(object, 
                    dims = NULL,
                    dimred.type = "PCA",
                    return.model = FALSE, 
                    umap.method = "umap", 
                    dimred.name = "UMAP", 
                    ...) {
             standardGeneric("RunUMAP")
           })

#' @export
setGeneric("RunTSNE", signature = "object",
           function(object,
                    dims = NULL, 
                    dimred.type = "PCA",
                    perplexity = 30, 
                    dimred.name = "TSNE", 
                    ...) {
             standardGeneric("RunTSNE")
           })

#' @export
setGeneric("GeneScatterPlot", signature = "object",
           function(object,
                    genes = "",
                    return.plot = FALSE,
                    dim.reduction.type = "tsne",
                    point.size = 0.7,
                    title = "",
                    plot.expressing.cells.last = FALSE,
                    nrow = NULL,
                    ncol = NULL) {
             standardGeneric("GeneScatterPlot")
           })

#' @export
setGeneric("FindAllClusterMarkers", signature = "object",
           function(object,
                    clustering.label,
                    test = "wilcox",
                    log2fc.threshold = 0.25,
                    min.pct = 0.1,
                    min.diff.pct = NULL,
                    min.cells.group = 3,
                    max.cells.per.cluster = NULL,
                    return.thresh = 0.01,
                    only.pos = FALSE) {
             standardGeneric("FindAllClusterMarkers")
           })

#' @export
setGeneric("FindClusterMarkers", signature = "object",
           function(object,
                    clustering.label,
                    clusters.1 = NULL,
                    clusters.2 = NULL,
                    test = "wilcox",
                    log2fc.threshold = 0.25,
                    min.pct = 0.1,
                    min.diff.pct = NULL,
                    min.cells.group = 3,
                    max.cells.per.cluster = NULL,
                    return.thresh = 0.01,
                    only.pos = FALSE) {
             standardGeneric("FindClusterMarkers")
           })

#' @export
setGeneric("VlnPlot", signature = "object",
           function(object,
                    clustering.label,
                    features,
                    return.plot = FALSE,
                    rotate.x.axis.labels = FALSE) {
             standardGeneric("VlnPlot")
           })

#' @export
setGeneric("HeatmapFeatures", signature = "object",
           function(object,
                    clustering.label,
                    features, 
                    ...) {
             standardGeneric("HeatmapFeatures")
           })

#' @export
setGeneric("GetCellClusterProbability", signature = "object",
           function(object,
                    icp.run = NULL,
                    icp.round = NULL,
                    concatenate = TRUE) {
               standardGeneric("GetCellClusterProbability")
           })

#' @export
setGeneric("SummariseCellClusterProbability", signature = "object",
           function(object,
                    icp.run = NULL,
                    icp.round = NULL,
                    funs = c("mean", "median"), 
                    scale.funs = TRUE, 
                    save.in.sce = TRUE) {
               standardGeneric("SummariseCellClusterProbability")
           })

#' @export
setGeneric("GetFeatureCoefficients", signature = "object",
           function(object,
                    icp.run = NULL,
                    icp.round = NULL) {
               standardGeneric("GetFeatureCoefficients")
           })

#' @export
setGeneric("MajorityVotingFeatures", signature = "object",
           function(object, label) {
               standardGeneric("MajorityVotingFeatures")
           })

#' @export
setGeneric("RunParallelDivisiveICP", signature = "object",
           function(object, batch.label = NULL, 
                    k = 16, d = 0.3, L = 50, 
                    r = 5, C = 0.3, reg.type = "L1", 
                    max.iter = 200, threads = 0,
                    icp.batch.size = Inf, 
                    train.with.bnn = TRUE, 
                    train.k.nn = 10, 
                    train.k.nn.prop = 0.3,
                    build.train.set = TRUE, 
                    build.train.params = list(),
                    scale.by = NULL, 
                    use.cluster.seed = TRUE,
                    divisive.method = "cluster.batch",
                    allow.free.k = TRUE,
                    ari.cutoff = 0.3,
                    verbose = FALSE) {
               standardGeneric("RunParallelDivisiveICP")
           })

#' @export
setGeneric("ReferenceMapping", signature = c("ref", "query"),
           function(ref, query, ref.label,
                    scale.query.by = NULL, 
                    project.umap = FALSE, 
                    select.icp.models = metadata(ref)$coralysis$pca.params$select.icp.tables, 
                    k.nn = 10, 
                    dimred.name.prefix = "") {
               standardGeneric("ReferenceMapping")
           })

#' @export
setGeneric("PlotDimRed", signature = "object",
           function(object, color.by, dimred = tail(reducedDimNames(object), n=1), 
                    dims = 1:2, use.color = NULL,  
                    point.size = 1, point.stroke = 1, legend.nrow = 2, seed.color = 123, 
                    label = FALSE, plot.theme = theme_classic(), 
                    rasterise = (ncol(object)<=3e4), rasterise.dpi = 300, 
                    legend.justification = "center", legend.size = 10, legend.title = color.by) {
               standardGeneric("PlotDimRed")
           })

#' @export
setGeneric("PlotExpression", signature = "object",
           function(object, color.by, dimred = tail(reducedDimNames(object), n=1), 
                    scale.values = FALSE, color.scale = "inferno", plot.theme = theme_classic(),
                    legend.title = color.by, point.size = 1, point.stroke = 1) {
               standardGeneric("PlotExpression")
           })

#' @export
setGeneric("PlotClusterTree", signature = "object",
           function(object, icp.run, color.by = NULL, use.color = NULL, 
                    seed.color = 123, legend.title = color.by, 
                    return.data = FALSE) {
               standardGeneric("PlotClusterTree")
           })

#' @export
setGeneric("BinCellClusterProbability", signature = "object",
           function(object, label, icp.run = NULL, icp.round = NULL, 
                    funs = "mean", bins = 20, aggregate.bins.by = "mean", 
                    use.assay = "logcounts") {
               standardGeneric("BinCellClusterProbability")
           })

#' @export
setGeneric("CellClusterProbabilityDistribution", signature = "object",
           function(object, label, group, probability = "scaled_mean_probs") {
               standardGeneric("CellClusterProbabilityDistribution")
           })

#' @export
setGeneric("TabulateCellBinsByGroup", signature = "object",
           function(object, group, relative = FALSE, margin = 1) {
               standardGeneric("TabulateCellBinsByGroup")
           })

#' @export
setGeneric("CellBinsFeatureCorrelation", signature = "object",
           function(object, labels = NULL, method = "pearson") {
               standardGeneric("CellBinsFeatureCorrelation")
           })