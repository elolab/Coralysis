#' @title Heatmap visualization of the gene markers identified by FindAllGeneMarkers
#'
#' @description
#' The GeneHeatmap function enables drawing a heatmap of the gene markers
#' identified by FindAllGeneMarkers, where the cell are grouped
#' by the clustering.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.type "manual" or "optimal". "manual" refers to the
#' clustering formed using the "SelectKClusters" function and "optimal"
#' to the clustering using the "CalcSilhInfo" function.
#' Default is "manual".
#' @param gene.markers a data frame of the gene markers generated
#' by FindAllGeneMarkers function. To accelerate the drawing, filtering
#' the dataframe by selecting e.g. top 10 genes is recommended.
#'
#' @name GeneHeatmap
#'
#' @return nothing
#'
#' @keywords gene heatmap grouped
#'
#' @importFrom S4Vectors metadata
#' @import pheatmap
#' @importFrom SingleCellExperiment logcounts
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareData(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,r=1,k=5) # Use L=200
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' gene_markers <- FindAllGeneMarkers(sce,log2fc.threshold = 0.5,min.pct = 0.5)
#' top10_log2FC <- SelectTopGenes(gene_markers,top.N=10,
#' criterion.type="log2FC",inverse=FALSE)
#' GeneHeatmap(sce,clustering.type = "manual",
#'  gene.markers = top10_log2FC)
#'
GeneHeatmap.SingleCellExperiment <- function(object,
                                             clustering.type,
                                             gene.markers) {
    
    
    if (clustering.type=="manual")
    {
        clustering <- metadata(object)$coralysis$clustering.manual
    } else if (clustering.type=="optimal")
    {
        clustering <- metadata(object)$coralysis$clustering.optimal
    } else {
        clustering <- metadata(object)$coralysis$clustering.manual
        cat("clustering.type='manual'")
    }
    
    data <- logcounts(object)
    data <- data[unique(gene.markers$gene),]
    # data <- scale(data,center = TRUE,scale = TRUE)
    data <- data[,order(clustering)]
    
    
    # Generate column annotations
    annotation = data.frame(cluster=sort(clustering))
    
    pheatmap(data,show_colnames = FALSE,
             gaps_col = cumsum(table(clustering[order(clustering)])),
             gaps_row = cumsum(table(gene.markers[!duplicated(gene.markers$gene),"cluster"])),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             annotation_col = annotation)
    
}
#' @rdname GeneHeatmap
#' @aliases GeneHeatmap
setMethod("GeneHeatmap", signature(object = "SingleCellExperiment"),
          GeneHeatmap.SingleCellExperiment)


#' @title Gene expression visualization using violin plots
#'
#' @description
#' The VlnPlot function enables visualizing expression levels of a gene,
#' or multiple genes, across clusters using Violin plots.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.type "manual" or "optimal". "manual"
#' refers to the clustering formed using the "SelectKClusters" function
#' and "optimal" to the clustering formed using the
#' "CalcSilhInfo" function. Default is "manual".
#' @param genes a character vector denoting the gene names that are visualized
#' @param return.plot return.plot whether to return the ggplot2 object
#' @param rotate.x.axis.labels a logical denoting whether the x-axis
#' labels should be rotated 90 degrees.
#' or just draw it. Default is \code{FALSE}.
#'
#' @name VlnPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords violin plot
#'
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom SingleCellExperiment logcounts
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareData(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' VlnPlot(sce,genes=c("CD3D","CD79A","CST3"))
#'
VlnPlot.SingleCellExperiment <- function(object,
                                         clustering.type,
                                         genes,
                                         return.plot,
                                         rotate.x.axis.labels) {
    
    
    if (clustering.type=="manual")
    {
        clustering <- metadata(object)$coralysis$clustering.manual
    } else if (clustering.type=="optimal")
    {
        clustering <- metadata(object)$coralysis$clustering.optimal
    } else {
        clustering <- metadata(object)$coralysis$clustering.manual
        message("clustering.type='manual'")
    }
    
    data <- logcounts(object)
    
    df <- as.numeric(t(data[genes,]))
    df <- data.frame(matrix(df,ncol = 1,dimnames = list(seq_len(length(df)),"Expression")))
    df$gene  <- unlist(lapply(genes,function(x) rep(x,ncol(data))))
    df$gene <- factor(df$gene)
    df$Cluster <- rep(as.character(clustering),length(genes))
    df$Cluster <- factor(df$Cluster)
    
    
    if (rotate.x.axis.labels)
    {
        plotlist <- lapply(genes,function(x) ggplot(df[df$gene==x,], aes_string(x='Cluster', y='Expression', fill='Cluster'))+geom_violin(trim=TRUE)+geom_jitter(height = 0, width = 0.1)+theme_classic()+ggtitle(x)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)))
    } else {
        plotlist <- lapply(genes,function(x) ggplot(df[df$gene==x,], aes_string(x='Cluster', y='Expression', fill='Cluster'))+geom_violin(trim=TRUE)+geom_jitter(height = 0, width = 0.1)+theme_classic()+ggtitle(x)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none"))
        
    }
    
    
    p <- plot_grid(plotlist = plotlist)
    
    if (return.plot)
    {
        return(p)
    } else {
        print(p)
    }
    
}
#' @rdname VlnPlot
#' @aliases VlnPlot
setMethod("VlnPlot", signature(object = "SingleCellExperiment"),
          VlnPlot.SingleCellExperiment)

#' @title Visualize gene expression over nonlinear dimensionality reduction
#'
#' @description
#' GeneScatterPlot enables visualizing gene expression of a gene over
#' nonlinear dimensionality reduction with t-SNE or UMAP.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param genes a character vector of the genes to be visualized
#' @param return.plot whether to return the ggplot2 object or just
#' draw it (default \code{FALSE})
#' @param dim.reduction.type "tsne" or "umap" (default "tsne")
#' @param point.size point size (default 0.7)
#' @param title text to write above the plot
#' @param plot.expressing.cells.last whether to plot the expressing genes
#' last to make the points more visible
#' @param nrow a positive integer that specifies the number of rows in
#' the plot grid. Default is \code{NULL}.
#' @param ncol a positive integer that specifies the number of columns
#' in the plot grid. Default is \code{NULL}.
#'
#' @name GeneScatterPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords gene scatter plot visualization
#'
#' @importFrom SingleCellExperiment reducedDim logcounts
#' @importFrom S4Vectors metadata metadata<-
#' @import ggplot2
#' @importFrom scales muted
#' @importFrom cowplot plot_grid
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareData(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- RunTSNE(sce)
#' GeneScatterPlot(sce,"CD14",dim.reduction.type="tsne")
#' sce <- RunUMAP(sce)
#' GeneScatterPlot(sce,"CD14",dim.reduction.type="umap")
#'
GeneScatterPlot.SingleCellExperiment <- function(object,
                                                 genes,
                                                 return.plot,
                                                 dim.reduction.type,
                                                 point.size,
                                                 title,
                                                 plot.expressing.cells.last,
                                                 nrow,
                                                 ncol) {
    
    if (dim.reduction.type=="umap") {
        two.dim.data <- reducedDim(object,"UMAP")
        xlab <- "UMAP_1"
        ylab <- "UMAP_2"
    } else if (dim.reduction.type=="tsne"){
        two.dim.data <- reducedDim(object,"TSNE")
        xlab <- "tSNE_1"
        ylab <- "tSNE_2"
    } else {
        stop("dim.reduction.type must be either 'tsne' or 'umap'")
    }
    
    if (length(genes)==1)
    {
        
        df <- as.data.frame(two.dim.data)
        
        if (!(genes %in% rownames(object))) {
            stop("invalid gene name")
        }
        
        color.by <- logcounts(object)[genes,]
        df$group <- color.by
        colnames(df) <- c("dim1","dim2","group")
        
        # Plot
        if (plot.expressing.cells.last) {
            df <- df[order(df$group,decreasing = FALSE),]
        }
        p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
            geom_point(size=point.size,aes_string(color='group')) +
            scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                   high = "blue",name = genes) +
            xlab(xlab) +
            ylab(ylab) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))
        
        if (title != "") {
            p <- p + ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5))
        }
        if (return.plot) {
            return(p)
        } else {
            print(p)
        }
    } else {
        
        plot_list <- list()
        for (gene in genes)
        {
            df <- as.data.frame(two.dim.data)
            
            if (!(gene %in% rownames(object)))
            {
                stop(paste0("invalid gene name: ",gene))
            }
            
            color.by <- logcounts(object)[gene,]
            df$group <- color.by
            colnames(df) <- c("dim1","dim2","group")
            
            # Plot 
            if (plot.expressing.cells.last) {
                df <- df[order(df$group,decreasing = FALSE),]
            }
            p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
                geom_point(size=point.size,aes_string(color='group')) +
                scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                       high = "blue",name = gene) +
                xlab(xlab) +
                ylab(ylab) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
            if (title != "") {
                p <- p + ggtitle(title) +
                    theme(plot.title = element_text(hjust = 0.5))
            }
            plot_list[[gene]] <- p
        }
        
        p <- plot_grid(plotlist = plot_list,align = "hv",nrow = nrow, ncol = ncol)
        
        if (return.plot) {
            return(p)
        } else {
            print(p)
        }
    }
    
}
#' @rdname GeneScatterPlot
#' @aliases GeneScatterPlot
setMethod("GeneScatterPlot", signature(object = "SingleCellExperiment"),
          GeneScatterPlot.SingleCellExperiment)


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
