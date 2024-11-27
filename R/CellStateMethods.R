#' @title Bin cell cluster probability
#' 
#' @description Bin cell cluster probability by a given cell label.
#' 
#' @param object An object of \code{SingleCellExperiment} class with ICP cell 
#' cluster probability tables saved in \code{metadata(object)$iloreg$joint.probability}. 
#' After running one of \code{RunParallelICP} or \code{RunParallelDivisiveICP}. 
#' @param label . 
#' @param icp.run ICP run(s) to retrieve from \code{metadata(object)$iloreg$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Specify a numeric vector to 
#' retrieve a specific set of tables. 
#' @param icp.round ICP round(s) to retrieve from \code{metadata(object)$iloreg$joint.probability}. 
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
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom dplyr %>% group_by mutate ntile ungroup
#'
BinCellClusterProbability.SingleCellExperiment <- function(object, label, icp.run, icp.round, 
                                                           funs, bins, aggregate.bins.by, 
                                                           use.assay) {
    # Retrieve important params
    L <- metadata(object)$iloreg$L
    k <- metadata(object)$iloreg$k
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), is(metadata(object)$iloreg$joint.probability, "list"), 
              is.numeric(L), is.numeric(k), any(is.null(icp.run), (is.numeric(icp.run) && all(icp.run <= L))), 
              any(is.null(icp.round), (is.numeric(icp.round))), 
              ((length(funs)==1) && is.character(funs) && (funs %in% c("mean", "median"))), 
              ((length(aggregate.bins.by)==1) && is.character(aggregate.bins.by) && (aggregate.bins.by %in% c("mean", "median"))), 
              use.assay %in% assayNames(object))
    
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
        mutate("probability_bins" = ntile(probability, n = bins)) %>% 
        mutate("probability_bins" = as.factor(probability_bins)) %>%
        group_by(label, probability_bins) %>% 
        mutate("aggregated_probability_bins" = get(aggregate.bins.by)(probability), 
               "aggregated_label_bins" = paste(label, paste0("bin", probability_bins), sep = "_")) %>% 
        ungroup(.)
    # Gene expression averaged by 'label' x bins
    gexp.bins <- AggregateClusterExpression(mtx = assay(x = object, i = use.assay), 
                                            cluster = bins.by.label$aggregated_label_bins, 
                                            fun = "mean")
    colnames(gexp.bins) <- gsub(pattern = "^cluster", replacement = "", x = colnames(gexp.bins))
    
    # SCE object
    col.data <- bins.by.label[,c("label", "probability_bins", "aggregated_probability_bins", "aggregated_label_bins")] %>% 
        distinct(.) %>% 
        as.data.frame(.) %>% 
        `row.names<-`(.$aggregated_label_bins)
    col.data <- col.data[colnames(gexp.bins),]
    sce <- SingleCellExperiment(assays = list("exp" = gexp.bins), 
                                colData = DataFrame(col.data), 
                                metadata = list("iloreg" = bins.by.label))
    assayNames(sce) <- use.assay
    return(sce)
}
#' @rdname BinCellClusterProbability
#' @aliases BinCellClusterProbability
setMethod("BinCellClusterProbability", signature(object = "SingleCellExperiment"),
          BinCellClusterProbability.SingleCellExperiment)
