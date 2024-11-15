#' @export
setGeneric("PrepareILoReg2",signature = "object",
           function(object) {
             standardGeneric("PrepareILoReg2")
           })

#' @export
setGeneric("RunParallelICP", signature = "object",
           function(object, batch.label = NULL, 
                    k = 15, d = 0.3, L = 200, 
                    r = 5, C = 0.3, reg.type = "L1", 
                    max.iter = 200, threads = 0,
                    icp.batch.size = Inf, 
                    train.with.bnn = TRUE, 
                    train.k.nn = 10, 
                    train.k.nn.prop = NULL,
                    build.train.set = TRUE, 
                    build.train.params = list(),
                    scale = FALSE, 
                    verbose = FALSE) {
             standardGeneric("RunParallelICP")
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
                    scale = FALSE,
                    center = TRUE,
                    threshold = 0,
                    method = "RSpectra", 
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
setGeneric("HierarchicalClustering", signature = "object",
           function(object) {
             standardGeneric("HierarchicalClustering")
           })

#' @export
setGeneric("CalcSilhInfo", signature = "object",
           function(object,
                    K.start = 2, K.end = 50) {
             standardGeneric("CalcSilhInfo")
           })

#' @export
setGeneric("SilhouetteCurve", signature = "object",
           function(object,
                    return.plot = FALSE) {
             standardGeneric("SilhouetteCurve")
           })

#' @export
setGeneric("SelectKClusters", signature = "object",
           function(object,
                    K = NULL) {
             standardGeneric("SelectKClusters")
           })

#' @export
setGeneric("MergeClusters", signature = "object",
           function(object,
                    clusters.to.merge = "",
                    new.name = "") {
             standardGeneric("MergeClusters")
           })

#' @export
setGeneric("RenameAllClusters", signature = "object",
           function(object,
                    new.cluster.names = "") {
             standardGeneric("RenameAllClusters")
           })

#' @export
setGeneric("RenameCluster", signature = "object",
           function(object,
                    old.cluster.name = "",
                    new.cluster.name = "") {
             standardGeneric("RenameCluster")
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
setGeneric("ClusteringScatterPlot",signature = "object",
           function(object,
                    clustering.type = "manual",
                    return.plot = FALSE,
                    dim.reduction.type = "",
                    point.size = 0.7,
                    title = "",
                    show.legend = TRUE) {
             standardGeneric("ClusteringScatterPlot")
           })

#' @export
setGeneric("FindAllGeneMarkers", signature = "object",
           function(object,
                    clustering.type = "manual",
                    test = "wilcox",
                    log2fc.threshold = 0.25,
                    min.pct = 0.1,
                    min.diff.pct = NULL,
                    min.cells.group = 3,
                    max.cells.per.cluster = NULL,
                    return.thresh = 0.01,
                    only.pos = FALSE) {
             standardGeneric("FindAllGeneMarkers")
           })

#' @export
setGeneric("FindGeneMarkers", signature = "object",
           function(object,
                    clusters.1 = NULL,
                    clusters.2 = NULL,
                    clustering.type = "",
                    test = "wilcox",
                    logfc.threshold = 0.25,
                    min.pct = 0.1,
                    min.diff.pct = NULL,
                    min.cells.group = 3,
                    max.cells.per.cluster = NULL,
                    return.thresh = 0.01,
                    only.pos = FALSE) {
             standardGeneric("FindGeneMarkers")
           })

#' @export
setGeneric("VlnPlot", signature = "object",
           function(object,
                    clustering.type = "manual",
                    genes = NULL,
                    return.plot = FALSE,
                    rotate.x.axis.labels=FALSE) {
             standardGeneric("VlnPlot")
           })

#' @export
setGeneric("GeneHeatmap", signature = "object",
           function(object,
                    clustering.type = "manual",
                    gene.markers = NULL) {
             standardGeneric("GeneHeatmap")
           })

#' @export
setGeneric("AnnotationScatterPlot", signature = "object",
           function(object,
                    annotation = NULL,
                    return.plot = FALSE,
                    dim.reduction.type = "",
                    point.size = 0.7,
                    show.legend = FALSE) {
             standardGeneric("AnnotationScatterPlot")
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
                    select.icp.models = metadata(ref)$iloreg$pca.params$select.icp.tables, 
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