#' @title Bin cell cluster probability
#' 
#' @description Bin cell cluster probability by a given cell label.
#' 
#' @param object An object of \code{SingleCellExperiment} class with ICP cell 
#' cluster probability tables saved in \code{metadata(object)$coralysis$joint.probability}. 
#' After running one of \code{RunParallelICP} or \code{RunParallelDivisiveICP}. 
#' @param label Label of interest available in \code{colData(object)} to group by the 
#' bins of cell cluster probability. 
#' @param icp.run ICP run(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Specify a numeric vector to 
#' retrieve a specific set of tables. 
#' @param icp.round ICP round(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Only relevant if probabilities
#' were obtained with the function \code{RunParallelDivisiveICP}, i.e., divisive ICP
#' was performed. Otherwise it is ignored and internally assumed as \code{icp.round = 1}, 
#' i.e., only one round. 
#' @param funs One function to summarise ICP cell cluster probability. One of \code{"mean"} 
#' or \code{"median"}. By default \code{"mean"}. 
#' @param bins Number of bins to bin cell cluster probability by cell \code{label} given. 
#' By default \code{20}. 
#' @param aggregate.bins.by One function to aggregate One of \code{"mean"} or 
#' \code{"median"}. By default \code{"mean"}.
#' @param use.assay Name of the assay that should be used to obtain the average expression 
#' of features across cell \code{label} probability bins. 
#' 
#' @name BinCellClusterProbability
#' 
#' @return A \code{SingleCellExperiment} class object with feature average expression by 
#' cell \code{label} probability bins. 
#' 
#' @keywords Bin cell cluster probability
#'
#' @importFrom SingleCellExperiment SingleCellExperiment 
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData assay assayNames assayNames<-
#' @importFrom dplyr %>% group_by mutate ntile ungroup distinct
#'
#' @examples 
#' \donttest{
#' # Packages
#' suppressPackageStartupMessages(library("SingleCellExperiment"))
#' 
#' # Import data from Zenodo
#' data.url <- "https://zenodo.org/records/14845751/files/pbmc_10Xassays.rds?download=1"
#' sce <- readRDS(file = url(data.url))
#' 
#' # Prepare data
#' sce <- PrepareData(object = sce)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "batch", L = 4, 
#'                               threads = 2)
#' 
#' # Cell states SCE object for a given cell type annotation or clustering
#' cellstate.sce <- BinCellClusterProbability(object = sce, label = "cell_type", 
#'                                            icp.round = 4, bins = 20) 
#' cellstate.sce
#' }
#'
BinCellClusterProbability.SingleCellExperiment <- function(object, label, icp.run, icp.round, 
                                                           funs, bins, aggregate.bins.by, 
                                                           use.assay) {
    # Retrieve important params
    L <- metadata(object)$coralysis$L
    k <- metadata(object)$coralysis$k
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), is(metadata(object)$coralysis$joint.probability, "list"), 
              is.numeric(L), is.numeric(k), any(is.null(icp.run), (is.numeric(icp.run) && all(icp.run <= L))), 
              any(is.null(icp.round), (is.numeric(icp.round))), 
              ((length(funs)==1) && is.character(funs) && (funs %in% c("mean", "median"))), 
              ((length(aggregate.bins.by)==1) && is.character(aggregate.bins.by) && (aggregate.bins.by %in% c("mean", "median"))), 
              use.assay %in% assayNames(object), 
              (is.character(label) && (length(label)==1) && (label %in% colnames(colData(object)))))
    
    # Getting probability
    probs <- SummariseCellClusterProbability(object = object, icp.run = icp.run, icp.round = icp.round, 
                                             funs = funs, scale.funs = FALSE, save.in.sce = FALSE)
    # Build data frame 
    data.bins <- data.frame(
        "cell_id" = colnames(object), 
        "label" = colData(object)[, label, drop=TRUE], 
        "probability" = probs[,paste0(funs, "_probs")] 
    )
    # Obtain probability bins by 'label'
    bins.by.label <- data.bins %>% 
        group_by(label) %>% 
        mutate("probability_bins" = ntile(.data$probability, n = bins)) %>% 
        mutate("probability_bins" = as.factor(.data$probability_bins)) %>%
        group_by(.data$label, .data$probability_bins) %>% 
        mutate("aggregated_probability_bins" = get(aggregate.bins.by)(.data$probability), 
               "aggregated_label_bins" = paste(.data$label, paste0("bin", .data$probability_bins), sep = "_")) %>% 
        ungroup()
    # Feature expression averaged by 'label' x bins
    gexp.bins <- AggregateClusterExpression(mtx = assay(x = object, i = use.assay), 
                                            cluster = bins.by.label$aggregated_label_bins, 
                                            fun = "mean")
    colnames(gexp.bins) <- gsub(pattern = "^cluster", replacement = "", x = colnames(gexp.bins))
    
    # SCE object
    col.data <- bins.by.label[,c("label", "probability_bins", "aggregated_probability_bins", "aggregated_label_bins")] %>% 
        distinct() %>% 
        as.data.frame() 
    row.names(col.data) <- col.data$aggregated_label_bins
    col.data <- col.data[colnames(gexp.bins),]
    bins.by.label <- as.data.frame(cbind(colData(object), bins.by.label))
    sce <- SingleCellExperiment(assays = list("exp" = gexp.bins), 
                                colData = DataFrame(col.data), 
                                metadata = list("coralysis" = bins.by.label))
    assayNames(sce) <- use.assay
    return(sce)
}
#' @rdname BinCellClusterProbability
#' @aliases BinCellClusterProbability
setMethod("BinCellClusterProbability", signature(object = "SingleCellExperiment"),
          BinCellClusterProbability.SingleCellExperiment)


#' @title Cell cluster probability distribution
#' 
#' @description Plot cell cluster probability distribution per label by group.
#' 
#' @param object An object of \code{SingleCellExperiment} class with aggregated 
#' cell cluster probability available in \code{colData(object)}, which can be 
#' obtained after running \code{SummariseCellClusterProbability()}. 
#' @param label Character specifying the \code{colData} variable to use as cell 
#' type/cluster label. 
#' @param group Character specifying the \code{colData} variable to use as 
#' categorical group variable.
#' @param probability Character specifying the aggregated cell cluster probability 
#' variable available in \code{colData}, used to plot its distribution. One of 
#' \code{"mean_probs"}, \code{"scaled_mean_probs"}, \code{"median_probs"}, 
#' \code{"scaled_median_probs"}. The availability of these variables in \code{colData} 
#' depends on the parameters given to the function \code{SummariseCellClusterProbability()} 
#' beforehand. By default assumes that \code{"scaled_mean_probs"} is available in
#' \code{colData}, which is only true if \code{SummariseCellClusterProbability()} 
#' function was run with \code{funs = "mean"} and \code{scale.funs = TRUE}. 
#' 
#' @name CellClusterProbabilityDistribution
#' 
#' @return A plot of class \code{ggplot}.
#' 
#' @keywords Distribution cell cluster probability
#'
#' @import ggplot2
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#'
#' @examples 
#' # Import package
#' suppressPackageStartupMessages(library("SingleCellExperiment"))
#' 
#' # Create toy SCE data
#' batches <- c("b1", "b2")
#' set.seed(239)
#' batch <- sample(x = batches, size = nrow(iris), replace = TRUE)
#' sce <- SingleCellExperiment(assays = list(logcounts = t(iris[,1:4])),  
#'                             colData = DataFrame("Species" = iris$Species, 
#'                                                 "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Prepare SCE object for analysis
#' sce <- PrepareData(sce)
#' 
#' # Multi-level integration (just for highlighting purposes; use default parameters)
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "Batch", 
#'                               k = 4, L = 25, C = 1, d = 0.5, 
#'                               train.with.bnn = FALSE, 
#'                               use.cluster.seed = FALSE,
#'                               build.train.set = FALSE, ari.cutoff = 0.1, 
#'                               threads = 2)
#' 
#' # Summarise cell cluster probability
#' sce <- SummariseCellClusterProbability(object = sce, icp.round = 2) # saved in 'colData'
#' 
#' # Search for differences in probabilities across group(s) 
#' # give an interesting variable to the "group" parameter
#' prob.dist <- CellClusterProbabilityDistribution(object = sce, label = "Species", 
#'                                                 group = "Batch", 
#'                                                 probability = "scaled_mean_probs")
#' prob.dist # print plot
#'
CellClusterProbabilityDistribution.SingleCellExperiment <- function(object, label, group, probability) {
    
    # Check input 
    stopifnot(is(object, "SingleCellExperiment"), is(metadata(object)$coralysis$joint.probability, "list"), 
              (is.character(label) && (length(label)==1) && (label %in% colnames(colData(object)))), 
              (is.character(group) && (length(group)==1) && (label %in% colnames(colData(object)))), 
              (is.character(probability) && (length(probability)==1) && 
                   (probability %in% c("mean_probs", "scaled_mean_probs", "median_probs", "scaled_median_probs")) && 
                   (probability %in% colnames(colData(object)))))
    
    # Get data
    data.plot <- as.data.frame(colData(object))
    
    # Plot
    p <- ggplot(data = data.plot, mapping = aes(x = .data[[probability]], color = .data[[group]])) + 
        stat_density(geom = "line", position = "identity", linewidth = 0.75) +
        facet_grid( ~ .data[[label]]) +
        theme_classic() + 
        scale_x_continuous(breaks = c(0, 0.5, 1)) + 
        theme(strip.background = element_rect(color = NA), 
              axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 7), 
              legend.position = "bottom", 
              legend.title.position = "top", 
              legend.title = element_text(hjust = 0.5))
    return(p)
} 
#' @rdname CellClusterProbabilityDistribution
#' @aliases CellClusterProbabilityDistribution
setMethod("CellClusterProbabilityDistribution", signature(object = "SingleCellExperiment"),
          CellClusterProbabilityDistribution.SingleCellExperiment)


#' @title Tabulate cell bins by group
#' 
#' @description Frequency of cells per cell cluster probability bin by group for each label. 
#' The label has to be specified beforehand to the function \code{BinCellClusterProbability()}.
#' 
#' @param object An object of \code{SingleCellExperiment} class obtained with the 
#' function \code{BinCellClusterProbability()}. 
#' @param group Character specifying the \code{colData} variable from the \code{SingleCellExperiment} 
#' object provided to the function \code{BinCellClusterProbability()} to use as 
#' categorical group variable.
#' @param relative Logical specifying if relative proportions of cell bins per 
#' \code{group} should be returned. By default \code{FALSE}, i.e., absolute values
#' are returned.  
#' @param margin If \code{relative} is \code{TRUE}, proportions should be calculated 
#' by: rows (\code{1}, the default); columns (\code{2}); or overall (\code{NULL}).
#' 
#' @name TabulateCellBinsByGroup
#' 
#' @return A list of tables with the frequency of cells per bin of cell cluster probability by group for each label. 
#' 
#' @keywords Table cell bins group
#'
#' @importFrom S4Vectors metadata
#'
#' @examples 
#' \donttest{
#' # Packages
#' suppressPackageStartupMessages(library("SingleCellExperiment"))
#' 
#' # Import data from Zenodo
#' data.url <- "https://zenodo.org/records/14845751/files/pbmc_10Xassays.rds?download=1"
#' sce <- readRDS(file = url(data.url))
#' 
#' # Prepare data
#' sce <- PrepareData(object = sce)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "batch", L = 4, 
#'                               threads = 2)
#' 
#' # Cell states SCE object for a given cell type annotation or clustering
#' cellstate.sce <- BinCellClusterProbability(object = sce, label = "cell_type", 
#'                                            icp.round = 4, bins = 20) 
#' cellstate.sce
#' 
#' # Tabulate cell bins by group 
#' # give an interesting variable to the "group" parameter
#' cellbins.tables <- TabulateCellBinsByGroup(object = cellstate.sce, 
#'                                            group = "batch", 
#'                                            relative = TRUE, 
#'                                            margin = 1)
#' }
#'
TabulateCellBinsByGroup.SingleCellExperiment <- function(object, group, relative, margin) {
    # Check input
    stopifnot(is(object, "SingleCellExperiment"), is(metadata(object)$coralysis, "data.frame"), 
              (is.character(group) && (length(group)==1) && (group %in% colnames(metadata(object)$coralysis))))
    
    # Parse data 
    col.data <- metadata(object)$coralysis
    
    # Split by label 
    col.data.by.label <- split(x = col.data, f = col.data$label)
    
    # Tabulate 
    cellbins.by.group <- lapply(X = col.data.by.label, FUN = function(x) {
        table(x[[group]], x$probability_bins)
    })
    if (relative) {
        cellbins.by.group <- lapply(X = cellbins.by.group, FUN = function(x) prop.table(x, margin = margin))
    }
    return(cellbins.by.group)
}
#' @rdname TabulateCellBinsByGroup
#' @aliases TabulateCellBinsByGroup
setMethod("TabulateCellBinsByGroup", signature(object = "SingleCellExperiment"),
          TabulateCellBinsByGroup.SingleCellExperiment)


#' @title Cell bins feature correlation
#' 
#' @description Correlation between cell bins for the given labels and features.
#' 
#' @param object An object of \code{SingleCellExperiment} class obtained with the 
#' function \code{BinCellClusterProbability()}. 
#' @param labels Character of label(s) from the label provided to the function 
#' \code{BinCellClusterProbability()}. By default \code{NULL}, i.e., all labels
#' are used. 
#' @param method Character specifying the correlation method to use. One of 
#' \code{"pearson"}, \code{"kendall"} or \code{"spearman"}. By default \code{"pearson"}
#' is used. 
#' 
#' @name CellBinsFeatureCorrelation
#' 
#' @return A data frame with the correlation coefficient for each feature (rows) 
#' across labels (columns). 
#' 
#' @keywords Cell bins feature correlation
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment logcounts
#'
#' @examples 
#' \donttest{
#' # Packages
#' suppressPackageStartupMessages(library("SingleCellExperiment"))
#' 
#' # Import data from Zenodo
#' data.url <- "https://zenodo.org/records/14845751/files/pbmc_10Xassays.rds?download=1"
#' sce <- readRDS(file = url(data.url))
#' 
#' # Prepare data
#' sce <- PrepareData(object = sce)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "batch", L = 4, 
#'                               threads = 2)
#' 
#' # Cell states SCE object for a given cell type annotation or clustering
#' cellstate.sce <- BinCellClusterProbability(object = sce, label = "cell_type", 
#'                                            icp.round = 4, bins = 20) 
#' cellstate.sce
#' 
#' # Pearson correlated features with "Monocyte"
#' cor.features.mono <- CellBinsFeatureCorrelation(object = cellstate.sce, 
#'                                                 labels = "Monocyte")
#' }
#'
CellBinsFeatureCorrelation.SingleCellExperiment <- function(object, labels, method) {
    
    # Check input
    stopifnot(is(object, "SingleCellExperiment"),
              (is.null(labels) || (is.character(labels) && all(labels %in% object$label))), 
              (is.character(method) && any(method %in% c("pearson", "kendall", "spearman"))))
    
    # Retrieve 'labels' if not given
    if (!is.factor(object$label)) {
        object$label <- as.factor(object$label)
    }
    if (is.null(labels)) {
        labels <- levels(object$label)
    }
    names(labels) <- labels
    
    # Subset data to 'labels'
    data.labels <- colData(object)
    data.labels <- data.labels[ data.labels$label %in% labels, ]
    data.labels <- split(x = as.data.frame(data.labels), f = data.labels$label)
    
    # Get correlation
    cor.features <- data.frame(matrix(ncol = length(labels), nrow = nrow(object)), row.names = row.names(object))
    colnames(cor.features) <- labels
    for (label in labels) {
        # Get probability & label bins per label
        bins <- data.labels[[label]][ order(data.labels[[label]]$probability_bins) , ]
        prop.bins <- bins$aggregated_probability_bins
        label.bins <- bins$aggregated_label_bins
        # Expression per cell bins per label
        avg.feature.exp <- logcounts(object[,label.bins])
        avg.feature.exp <- avg.feature.exp[ rowSums(avg.feature.exp>0)>0, ] # remove non-expressed features
        cor.feature.cell.bins <- apply(X = avg.feature.exp, MARGIN = 1, FUN = function(x) {
            stats::cor(x = x, y = prop.bins, method = method)
        })
        cor.features[names(cor.feature.cell.bins), label] <- cor.feature.cell.bins
    }
    return(cor.features)
}
#' @rdname CellBinsFeatureCorrelation
#' @aliases CellBinsFeatureCorrelation
setMethod("CellBinsFeatureCorrelation", signature(object = "SingleCellExperiment"),
          CellBinsFeatureCorrelation.SingleCellExperiment)
