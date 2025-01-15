#' @title Heatmap visualization of the expression of features by clusters
#'
#' @description The \code{HeatmapFeatures} function draws a heatmap of features 
#' by cluster identity.  
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.label A variable name (of class \code{character}) available 
#' in the cell metadata \code{colData(object)} with the clustering labels 
#' (\code{character} or \code{factor}) to use.
#' @param features Feature names to plot by cluster (\code{character}) matching 
#' \code{row.names(object)}.
#' @param ... Parameters to pass to \code{pheatmap::pheatmap} function.
#'
#' @name HeatmapFeatures
#'
#' @return nothing
#'
#' @keywords feature heatmap grouped
#'
#' @importFrom S4Vectors metadata
#' @import pheatmap
#' @importFrom SingleCellExperiment logcounts
#' @examples
#' # Import package
#' library("SingleCellExperiment")
#' 
#' # Prepare SCE object for analysis
#' pbmc_10Xassays <- PrepareData(pbmc_10Xassays)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes; use 'L=50' or greater
#' set.seed(123)
#' pbmc_10Xassays <- RunParallelDivisiveICP(object = pbmc_10Xassays, batch.label = "batch", L = 4, threads = 1) 
#' 
#' # Run PCA 
#' pbmc_10Xassays <- RunPCA(pbmc_10Xassays, p = 10)
#' 
#' # Plot features by cell type (give a clustering instead)
#' genes <- c("CD8A", "CCR7", "GZMA", "MS4A1", "CD27")
#' HeatmapFeatures(object = pbmc_10Xassays, clustering.label = "cell_type", features = genes) # log-normalized
#' HeatmapFeatures(object = pbmc_10Xassays, clustering.label = "cell_type", features = genes, scale = "row") # scale genes
#' 
HeatmapFeatures.SingleCellExperiment <- function(object, clustering.label, features, ...) {
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), 
              all(is.character(clustering.label), length(clustering.label)==1, clustering.label %in% colnames(colData(object)), (is.character(object[[clustering.label]]) || is.factor(object[[clustering.label]]))), 
              all(is.character(features), all(features %in% row.names(object))))
    
    # Retrieve clustering as factor
    clustering <- object[[clustering.label]]
    if (is.character(clustering)) { # if character, coerce to factor
        clustering <- as.factor(clustering)
    }
    names(clustering) <- colnames(object)
    clusters <- levels(clustering)
    
    # Retrieve data
    data <- logcounts(object)
    data <- data[unique(features),]
    data <- data[,order(clustering)]
    
    # Generate column annotations
    annotation = data.frame(cluster=sort(clustering))
    
    # Plot heatmap
    extra.params <- list(...)
    extra.params$show_colnames <- ifelse(is.null(extra.params$show_colnames), FALSE, extra.params$show_colnames)
    extra.params$cluster_rows <- ifelse(is.null(extra.params$cluster_rows), FALSE, extra.params$cluster_rows)
    extra.params$cluster_cols <- ifelse(is.null(extra.params$cluster_cols), FALSE, extra.params$cluster_cols)
    if (is.null(extra.params$gaps_col)) extra.params$gaps_col <- cumsum(table(clustering[order(clustering)]))
    if (is.null(extra.params$annotation_col)) extra.params$annotation_col <- annotation
    params <- c(list(mat = data), extra.params)
    do.call(pheatmap, params)
}
#' @rdname HeatmapFeatures
#' @aliases HeatmapFeatures
setMethod("HeatmapFeatures", signature(object = "SingleCellExperiment"),
          HeatmapFeatures.SingleCellExperiment)


#' @title Visualization of feature expression using violin plots
#'
#' @description The \code{VlnPlot} function enables visualizing expression levels 
#' of feature(s), across clusters using violin plots.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.label A variable name (of class \code{character}) available 
#' in the cell metadata \code{colData(object)} with the clustering labels 
#' (\code{character} or \code{factor}) to use.
#' @param features Feature names to plot by cluster (\code{character}) matching 
#' \code{row.names(object)}.
#' @param return.plot return.plot whether to return the \code{ggplot2} object. 
#' Default is \code{FALSE}.
#' @param rotate.x.axis.labels a logical denoting whether the x-axis labels should 
#' be rotated 90 degrees or just draw it. Default is \code{FALSE}.
#'
#' @name VlnPlot
#'
#' @return A \code{ggplot2} object if \code{return.plot=TRUE}.
#'
#' @keywords violin plot
#'
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom SingleCellExperiment logcounts
#'
#' @examples
#' # Import package
#' library("SingleCellExperiment")
#' 
#' # Prepare SCE object for analysis
#' pbmc_10Xassays <- PrepareData(pbmc_10Xassays)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes; use 'L=50' or greater
#' set.seed(123)
#' pbmc_10Xassays <- RunParallelDivisiveICP(object = pbmc_10Xassays, batch.label = "batch", L = 4, threads = 1) 
#' 
#' # Run PCA 
#' pbmc_10Xassays <- RunPCA(pbmc_10Xassays, p = 10)
#' 
#' # Plot features by cell type (give a clustering instead)
#' genes <- c("CD8A", "CCR7", "GZMA", "MS4A1", "CD27")
#' VlnPlot(pbmc_10Xassays, clustering.label = "cell_type", features = genes, rotate.x.axis.labels = TRUE)
#' 
VlnPlot.SingleCellExperiment <- function(object, clustering.label, features, return.plot, rotate.x.axis.labels) {
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), 
              all(is.character(clustering.label), length(clustering.label)==1, clustering.label %in% colnames(colData(object)), (is.character(object[[clustering.label]]) || is.factor(object[[clustering.label]]))), 
              all(is.character(features), all(features %in% row.names(object))), 
              is.logical(return.plot), 
              is.logical(rotate.x.axis.labels))
    
    # Retrieve clustering as factor
    clustering <- object[[clustering.label]]
    if (is.character(clustering)) { # if character, coerce to factor
        clustering <- as.factor(clustering)
    }
    names(clustering) <- colnames(object)
    clusters <- levels(clustering)
    
    # Retrieve data
    data <- logcounts(object)
    df <- as.numeric(t(data[features,]))
    df <- data.frame(matrix(df, ncol = 1, dimnames = list(seq_len(length(df)), "Expression")))
    df$feature  <- unlist(lapply(features,function(x) rep(x, ncol(data))))
    df$feature <- factor(df$feature)
    df$Cluster <- rep(as.character(clustering), length(features))
    df$Cluster <- factor(df$Cluster)
    
    # Plot violin plot(s)
    if (rotate.x.axis.labels) {
        plotlist <- lapply(features,function(x) ggplot(df[df$feature==x,], aes_string(x='Cluster', y='Expression', fill='Cluster')) + 
                               geom_violin(trim=TRUE) + 
                               geom_jitter(height = 0, width = 0.1) + 
                               theme_classic() + 
                               ggtitle(x) + 
                               theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
    } else {
        plotlist <- lapply(features,function(x) ggplot(df[df$feature==x,], aes_string(x='Cluster', y='Expression', fill='Cluster')) + 
                               geom_violin(trim = TRUE) + 
                               geom_jitter(height = 0, width = 0.1) + 
                               theme_classic() + 
                               ggtitle(x) + 
                               theme(plot.title = element_text(hjust = 0.5), legend.position = "none"))
    }
    p <- plot_grid(plotlist = plotlist)
    
    if (return.plot) {
        return(p)
    } else {
        print(p)
    }
}
#' @rdname VlnPlot
#' @aliases VlnPlot
setMethod("VlnPlot", signature(object = "SingleCellExperiment"),
          VlnPlot.SingleCellExperiment)


#' @title Plot dimensional reduction categorical variables
#' 
#' @description Plot categorical variables in dimensional reduction. 
#' 
#' @param object An object of \code{SingleCellExperiment} class.
#' @param color.by Categorical variable available in \code{colData(object)} to 
#' plot. 
#' @param dimred Dimensional reduction available in \code{ReducedDimNames(object)}
#' to plot. By default the last dimensional reduction in the object is used. 
#' @param dims Dimensions from the dimensional reduction embedding to plot.
#' @param use.color Character specifying the colors. By default \code{NULL}, i.e., 
#' colors are randomly chosen based on the seed given at \code{seed.color}. 
#' @param point.size Size of points. By default \code{1}. 
#' @param point.stroke Size of stroke. By default \code{1}. 
#' @param legend.nrow Display legend items by this number of rows. By default \code{2}.
#' @param seed.color Seed to randomly select colors. By default \code{123}.
#' @param label Logical to add or not categorical labels to the centroid categories. 
#' By default \code{FALSE}, i.e., labels are not added. 
#' @param plot.theme Plot theme available in \code{ggplot2}. By default \code{theme_classic()}. 
#' @param rasterise Logical specifying if points should be rasterised or not. By 
#' default \code{TRUE}, if more than 3e4 cells, otherwise \code{FALSE}. 
#' @param rasterise.dpi In case \code{rasterise = TRUE}, DPI to use. By default 
#' \code{300}. 
#' @param legend.justification Legend justification. By default \code{"center"}. 
#' @param legend.size Legend size. By default \code{10}
#' @param legend.title Legend title. By default the same as given at \code{color.by}. 
#' 
#' @name PlotDimRed
#' 
#' @return A plot of class \code{ggplot}. 
#' 
#' @keywords Dimensional reduction visualization
#'
#' @import ggplot2
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom dplyr %>% mutate mutate_at group_by
#'
#' @examples
#' # Import package
#' library("SingleCellExperiment")
#' 
#' # Prepare SCE object for analysis
#' pbmc_10Xassays <- PrepareData(pbmc_10Xassays)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes; use 'L=50' or greater
#' set.seed(123)
#' pbmc_10Xassays <- RunParallelDivisiveICP(object = pbmc_10Xassays, batch.label = "batch", L = 4, threads = 1) 
#' 
#' # Run PCA 
#' pbmc_10Xassays <- RunPCA(pbmc_10Xassays, p = 10)
#' 
#' # Run UMAP
#' set.seed(123)
#' pbmc_10Xassays <- RunUMAP(pbmc_10Xassays, dimred.type = "PCA")
#' 
#' # Plot batch 
#' PlotDimRed(object = pbmc_10Xassays, color.by = "batch", dimred = "UMAP")
#' 
#' # Plot cell type annotations
#' PlotDimRed(object = pbmc_10Xassays, color.by = "cell_type", legend.nrow = 3, dimred = "UMAP", label = TRUE)
#' 
PlotDimRed.SingleCellExperiment <- function(object, color.by, dimred, dims, use.color,  
                                            point.size, point.stroke, legend.nrow, seed.color, 
                                            label, plot.theme, rasterise, rasterise.dpi, 
                                            legend.justification, legend.size, legend.title) {
    
    # Check input
    stopifnot(is(object, "SingleCellExperiment"), all(is.character(color.by), color.by %in% colnames(colData(object))))
    
    # Parse data
    data.plot <- as.data.frame(reducedDim(object, dimred)[,dims])
    axes <- paste0(dimred, dims)
    colnames(data.plot) <- axes
    data.plot <- cbind(data.plot, colData(object))
    if (is.null(use.color)) {
        color.palettes <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
        color.palette <- unlist(mapply(RColorBrewer::brewer.pal, color.palettes$maxcolors, rownames(color.palettes)))
        data.plot[,color.by] <- as.factor(data.plot[,color.by])
        ngroups <- nlevels(data.plot[,color.by, drop=TRUE])
        set.seed(seed.color)
        use.color <- sample(color.palette, ngroups)
    }
    if (label) {
        data.plot <- data.plot %>% as.data.frame() %>% group_by(.data[[color.by]]) %>% 
            mutate(x = mean(.data[[axes[1]]]), y = mean(.data[[axes[2]]]), label = .data[[color.by]]) %>% 
            mutate_at(vars(label), list(~replace(., duplicated(.), NA)))
    }
    p <- ggplot(data = data.plot, mapping = aes(x = .data[[axes[1]]], y = .data[[axes[2]]], color = .data[[color.by]])) + 
        scale_color_manual(values = use.color) + 
        plot.theme + 
        theme(legend.position = "bottom", 
              legend.box.margin = margin(0, 0, 0, 0), 
              legend.key.size = unit(0.1, 'mm'), 
              legend.justification = legend.justification, 
              legend.text = element_text(size = legend.size), 
              legend.title.position = "top", legend.title = element_text(hjust = 0.5)) +
        guides(color = guide_legend(title=legend.title, nrow = legend.nrow, bycol = TRUE, 
                                    override.aes = list(size=2.5)))
    if (label) {
        p <- p + ggrepel::geom_text_repel(mapping = aes(x = x, y = y, label = label), 
                                          size = 3.5, nudge_x = .15, box.padding = 0.5,
                                          nudge_y = 1, segment.curvature = -0.1,
                                          segment.ncp = 3, segment.angle = 20, 
                                          segment.size = 0.05, 
                                          fontface = "bold")
    }
    if (rasterise) { # rasterise dots if >=30K cells
        p + ggrastr::rasterise(geom_point(size = point.size,  stroke = point.stroke), dpi = rasterise.dpi)
    } else {
        p + geom_point(size = point.size, stroke = point.stroke)
    }
}
#' @rdname PlotDimRed
#' @aliases PlotDimRed
setMethod("PlotDimRed", signature(object = "SingleCellExperiment"),
          PlotDimRed.SingleCellExperiment)


#' @title Plot dimensional reduction feature expression
#' 
#' @description Plot feature expression in dimensional reduction. 
#' 
#' @param object An object of \code{SingleCellExperiment} class.
#' @param color.by Categorical variable available in \code{colData(object)} to 
#' plot. 
#' @param dimred Dimensional reduction available in \code{ReducedDimNames(object)}
#' to plot. By default the last dimensional reduction in the object is used. 
#' @param color.scale Character of color scale palette to be passed to 
#' \code{ggplot2::scale_color_viridis_c}. By default \code{inferno}. Other palettes
#' are also available such as \code{viridis}.  
#' @param point.size Size of points. By default \code{1}. 
#' @param point.stroke Size of stroke. By default \code{1}. 
#' @param plot.theme Plot theme available in \code{ggplot2}. By default \code{theme_classic()}. 
#' @param legend.title Legend title. By default the same as given at \code{color.by}. 
#' 
#' @name PlotExpression
#' 
#' @return A plot of class \code{ggplot}. 
#' 
#' @keywords Dimensional reduction visualization
#'
#' @import ggplot2
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom dplyr %>% arrange
#'
#' @examples 
#' # Import package
#' library("SingleCellExperiment")
#' 
#' # Prepare SCE object for analysis
#' pbmc_10Xassays <- PrepareData(pbmc_10Xassays)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes; use 'L=50' or greater
#' set.seed(123)
#' pbmc_10Xassays <- RunParallelDivisiveICP(object = pbmc_10Xassays, batch.label = "batch", L = 4, threads = 1) 
#' 
#' # Run PCA 
#' pbmc_10Xassays <- RunPCA(pbmc_10Xassays, p = 10)
#' 
#' # Run UMAP
#' set.seed(123)
#' pbmc_10Xassays <- RunUMAP(pbmc_10Xassays, dimred.type = "PCA")
#' 
#' # Plot cell cluster probabilities - mean
#' # possible options: "mean_probs", "median_probs", "scaled_median_probs" 
#' pbmc_10Xassays <- SummariseCellClusterProbability(object = pbmc_10Xassays, icp.round = 4, scale.funs = TRUE, save.in.sce = TRUE) # save result in 'colData'
#' PlotExpression(object = pbmc_10Xassays, color.by = "scaled_mean_probs", dimred = "UMAP", color.scale = "viridis", legend.title = "Scaled mean prob.\n(min-max)") 
#' 
#' # Plot expression level of one or more features
#' ## one 
#' PlotExpression(object = pbmc_10Xassays, color.by = "LGALS2", dimred = "UMAP")
#' 
#' ## more than one
#' genes <- c("CD8A", "CCR7", "GZMA", "MS4A1", "CD27")
#' exp.plots <- lapply(X = genes, FUN = function(x) {
#'     PlotExpression(object = pbmc_10Xassays, color.by = x, dimred = "UMAP", scale.values = TRUE, point.size = 0.25, point.stroke = 0.25)
#' })
#' cowplot::plot_grid(plotlist = exp.plots, ncol = 5, align = "vh")
#' 
PlotExpression.SingleCellExperiment <- function(object, color.by, dimred, scale.values, color.scale, 
                                                plot.theme, legend.title, point.size, point.stroke) {
    
    # Check input
    stopifnot(is(object, "SingleCellExperiment"), all(is.character(color.by), color.by %in% colnames(colData(object)) || color.by %in% row.names(object)))
    
    # Parse data
    data.plt <- as.data.frame(colData(object))
    if (color.by %in% row.names(object)) {
        values <- logcounts(object)[color.by,]
        if (scale.values) {
            values <- scale(values)[,1]
        }
        data.plt <- cbind(data.plt, "color_by" = values)
        colnames(data.plt)[colnames(data.plt) == "color_by"] <- color.by 
    }
    dim.red <- reducedDim(object, type = dimred)[,1:2]
    axes <- paste0(dimred, 1:2)
    colnames(dim.red) <- axes
    data.plt <- cbind(data.plt, dim.red) %>% 
        arrange(.data[[color.by]])
    
    # Plot
    p <- ggplot(data = data.plt, mapping = aes(x = .data[[axes[1]]], y = .data[[axes[2]]], color = .data[[color.by]])) + 
        geom_point(size = point.size, stroke = point.stroke) + 
        plot.theme + 
        scale_color_viridis_c(name = ifelse(!is.null(legend.title), legend.title, color.by), option = color.scale)
    return(p)
}
#' @rdname PlotExpression
#' @aliases PlotExpression
setMethod("PlotExpression", signature(object = "SingleCellExperiment"),
          PlotExpression.SingleCellExperiment)


#' @title Plot cluster tree 
#' 
#' @description Plot cluster tree by or cluster probability or categorical variable.
#' 
#' @param object An object of \code{SingleCellExperiment} class.
#' @param icp.run ICP run(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Specify a numeric vector to 
#' retrieve a specific set of tables. 
#' @param color.by Categorical variable available in \code{colData(object)} to 
#' plot. If \code{NULL} the cluster probability is represented instead. 
#' By default \code{NULL}.  
#' @param use.color Character specifying the colors. By default \code{NULL}, i.e., 
#' colors are randomly chosen based on the seed given at \code{seed.color}. 
#' @param seed.color Seed to randomly select colors. By default \code{123}.
#' @param legend.title Legend title. By default the same as given at \code{color.by}.
#' Ignored if \code{color.by} is \code{NULL}.
#' @param return.data Return data frame used to plot. Logical. By default \code{FALSE}, 
#' i.e., only the plot is returned.
#' 
#' @name PlotClusterTree
#' 
#' @return A plot of class \code{ggplot} or a list with a plot of class \code{ggplot} and a data frame. 
#' 
#' @keywords Dimensional reduction visualization
#'
#' @import ggplot2
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom dplyr %>% mutate summarise group_by across all_of ungroup n
#'
#' @examples 
#' # Import package
#' library("SingleCellExperiment")
#' 
#' # Prepare SCE object for analysis
#' pbmc_10Xassays <- PrepareData(pbmc_10Xassays)
#' 
#' # Multi-level integration - 'L = 4' just for highlighting purposes; use 'L=50' or greater
#' set.seed(123)
#' pbmc_10Xassays <- RunParallelDivisiveICP(object = pbmc_10Xassays, batch.label = "batch", L = 4, threads = 1) 
#' 
#' # Plot probability
#' PlotClusterTree(object = pbmc_10Xassays, icp.run = 2)
#' 
#' # Plot batch label distribution
#' PlotClusterTree(object = pbmc_10Xassays, icp.run = 2, color.by = "batch")
#' 
#' # Plot cell type label distribution
#' PlotClusterTree(object = pbmc_10Xassays, icp.run = 2, color.by = "cell_type")
#' 
PlotClusterTree.SingleCellExperiment <- function(object, icp.run, color.by, use.color, 
                                                 seed.color, legend.title, return.data) {
    
    # Get probabilities for icp.run
    divisive.icp <- metadata(object)$coralysis$divisive.icp
    stopifnot(divisive.icp)
    k <- metadata(object)$coralysis$k
    icp.round <- 1:log2(k)
    probs <- GetCellClusterProbability(object = object, icp.run = icp.run, icp.round = icp.round, concatenate = FALSE)
    
    # Parse data to plot
    k.rounds.inverted <- rev(2**icp.round)
    node.pos <- node.seg <- list()
    for (i in seq_along(k.rounds.inverted)) {
        if (i == 1) {
            node.pos[[i]] <- 1:k.rounds.inverted[i]
        } else {
            node.pos[[i]] <- node.pos[[ii]][c(T, F)] + value
        }
        ii <- i
        value <- diff(node.pos[[ii]][1:2])/2
        node.seg[[i]] <- sort(c(node.pos[[ii]][c(T, F)] + value, node.pos[[ii]][c(F, T)] - value)) 
    }
    node.pos <- rev(node.pos)
    node.seg <- rev(node.seg)
    names(node.seg) <- names(node.pos) <- paste0("K", rev(k.rounds.inverted))
    
    # Get clusters & cell cluster probability
    cell.clusters <- SummariseCellClusterProbability(object = object, icp.run = icp.run, 
                                                     icp.round = icp.round, funs = NULL, 
                                                     save.in.sce = FALSE)
    
    # Parse data
    df.list <- list()
    col.names <- colnames(cell.clusters)
    for (i in icp.round) {
        df.list[[i]] <- data.frame(
            "K" = as.factor(i), 
            "k" = as.numeric(cell.clusters[,col.names[c(T,F)][i]]), 
            "p" = cell.clusters[,col.names[c(F,T)][i]]
        )
        df.list[[i]][,"start"] <- node.pos[[i]][df.list[[i]][,"k"]]
        df.list[[i]][,"end"] <- node.seg[[i]][df.list[[i]][,"k"]]
        df.list[[i]] <- cbind(df.list[[i]], colData(object))
    }
    df <- do.call(rbind, df.list)
    
    # Plot cluster tree
    if (is.null(color.by)) { # if 'color.by' is NULL, plot median cluster probability tree
        group.by <- c("K", "k", "start", "end")
        df <- df %>%
            group_by(., across(all_of(group.by))) %>% 
            summarise(., median = median(p), n = n(), .groups = "keep") %>% 
            ungroup(.) %>% 
            mutate(., Kannot = factor(paste0("K", 2**(icp.round))[K], levels = rev(paste0("K", 2**(icp.round)))))
        p <- ggplot(data = df, mapping = aes(x = start, y = Kannot)) +
            geom_segment(mapping = aes(x = start, xend = start, 
                                       yend = as.numeric(Kannot)+1, y = Kannot),
                         color = "black") +
            geom_segment(mapping = aes(x = start, xend = end, 
                                       y = as.numeric(Kannot)+1, yend = as.numeric(Kannot)+1), 
                         color = "black") + 
            geom_point(mapping = aes(size = n, color = median)) + 
            theme_minimal() + 
            theme(panel.grid.minor.x = element_blank(), 
                  panel.grid.major.x = element_blank(), 
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank(), 
                  legend.position = "bottom", 
                  legend.box = "horizontal", 
                  legend.title.position = "top", 
                  axis.text.y = element_text(size = 14, face = "italic"), 
                  plot.title = element_text(hjust = 0.5), 
                  legend.title.align = 0.5) + 
            scale_color_viridis_c(name = "Median Cluster Probability", limits = c(0, 1)) + 
            scale_size(name = "Cell Cluster Size") + 
            ggtitle(paste0("Divisive ICP Cluster Tree (L=", icp.run, ")"))
    } else { # plot cluster tree by 'color.by'
        group.by <- c("K", "k", "start", "end", color.by) 
        df <- df %>%
            group_by(., across(all_of(group.by))) %>% 
            summarise(., n = n(), .groups = "keep") %>% 
            ungroup(.) %>% 
            mutate(., Kannot = factor(paste0("K", 2**(icp.round))[K], levels = rev(paste0("K", 2**(icp.round))))) %>% 
            tidyr::pivot_wider(., names_from = .data[[color.by]], values_from = n)
        df[is.na(df)] <- 0
        object[[color.by]] <- as.factor(object[[color.by]])
        if (is.null(use.color)) {
            color.palettes <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
            color.palette <- unlist(mapply(RColorBrewer::brewer.pal, color.palettes$maxcolors, rownames(color.palettes)))
            ngroups <- nlevels(object[[color.by]])
            set.seed(seed.color)
            use.color <- sample(color.palette, ngroups)
        }
        p <- ggplot(data = df) +
            geom_segment(mapping = aes(x = start, xend = start, 
                                       yend = as.numeric(Kannot)+1, y = Kannot),
                         color = "black") +
            geom_segment(mapping = aes(x = start, xend = end, 
                                       y = as.numeric(Kannot)+1, yend = as.numeric(Kannot)+1), 
                         color = "black") + 
            scatterpie::geom_scatterpie(data = df, mapping = aes(x = start, y = as.numeric(Kannot), group = k), 
                                        cols = levels(object[[color.by]])) + 
            coord_equal() + 
            theme_minimal() + 
            theme(panel.grid.minor.x = element_blank(), 
                  panel.grid.major.x = element_blank(), 
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank(), 
                  legend.title.position = "top", 
                  legend.title = element_text(hjust = 0.5, size = 12.5),
                  axis.text.y = element_text(size = 14, face = "italic"), 
                  plot.title = element_text(hjust = 0.5), 
                  legend.position = "bottom", 
                  legend.box = "horizontal", 
                  #legend.title.align = 0.5, 
                  legend.text = element_text(size=12),
                  legend.spacing.y = unit(2, "cm")) + 
            ggtitle(paste0("Divisive ICP Cluster Tree (L=", icp.run, ")")) + 
            scale_fill_manual(name = legend.title, values = use.color) + 
            guides(fill = guide_legend(override.aes = list(size=3.5)))
    }
    if (return.data) {
        out <- list(data = df, plot = p)
        return(out)
    } else {
        return(p)
    }
}
#' @rdname PlotClusterTree
#' @aliases PlotClusterTree
setMethod("PlotClusterTree", signature(object = "SingleCellExperiment"),
          PlotClusterTree.SingleCellExperiment)
