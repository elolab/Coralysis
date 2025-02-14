#' @title Prepare \code{SingleCellExperiment} object for analysis
#'
#' @description This function prepares the \code{SingleCellExperiment} object 
#' for analysis. The only required input is an object of class \code{SingleCellExperiment} 
#' with at least data in the \code{logcounts} slot.
#'
#' @param object An object of \code{SingleCellExperiment} class.
#'
#' @name PrepareData
#'
#' @return An object of \code{SingleCellExperiment} class.
#'
#' @keywords prepare clean normalized data
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom methods is
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))

# # Prepare SCE object for analysis
#' sce <- PrepareData(sce)
#'
PrepareData.SingleCellExperiment <- function(object) {
    
    # Check that there are data in `logcounts` slot
    if (!("logcounts" %in% assayNames(object))) {
        stop(paste("`Error: `logcounts` slot is missing from your ",
                   "SingleCellExperiment object. This can be any kind of ",
                   "normalized data matrix. Set it by executing ",
                   "logcounts(object) <- norm_data",sep = ""))
        return(object)
    }
    
    # Remove duplicate features from the data in `logcounts` slot
    if (sum(duplicated(rownames(object))) != 0) {
        features_before <- length(rownames(object))
        object <- object[!duplicated(rownames(object)), ]
        features_after <- length(rownames(object))
        message(paste("data in SingleCellExperiment object contained duplicate ",
                      " features. ", features_before - features_after,
                      "/", features_before, " were filtered out."))
    }
    
    # Convert the data in `logcounts` slot into object of `dgCMatrix` class.
    if (is(logcounts(object), "matrix")) {
        logcounts(object) <- Matrix::Matrix(logcounts(object), sparse = TRUE)
        message(paste("Converting object of `matrix` class into `dgCMatrix`.",
                      " Please note that Coralysis has been designed to work with ",
                      "sparse data, i.e. data with ",
                      "a high proportion of zero values! Dense data will likely " ,
                      "increase run time and memory usage drastically!",sep=""))
    }
    else if (is(logcounts(object), "data.frame")) {
        logcounts(object) <- Matrix::Matrix(as.matrix(logcounts(object)), sparse = TRUE)
        message(paste("Converting object of `data.frame` class into `dgCMatrix`.",
                      " Please note that Coralysis has been designed to work with ",
                      "sparse data, i.e. data with ",
                      "a high proportion of zero values!",sep = ""))
    }
    else if (is(logcounts(object), "dgCMatrix")) {
        message("Data in `logcounts` slot already of `dgCMatrix` class...")
    }
    else {
        stop("Error: Data in `logcounts` slot is not of `matrix`, `data.frame` ",
             "or `dgCMatrix` class.")
        return(object)
    }
    
    # Filter features that are not expressed in any of the cells
    features_before_filtering <- nrow(object)
    non_expressing_features <- rownames(object)[Matrix::rowSums(logcounts(object)) != 0]
    object <- object[non_expressing_features,]
    features_after_filtering <- nrow(object)
    message(paste(features_after_filtering,"/",features_before_filtering,
                  " features remain after filtering features with only zero values.",
                  sep = ""))
    
    # Create a place into `metadata`` slot for the data from Coralysis
    metadata(object)$coralysis <- list()
    
    return(object)
}
#' @rdname PrepareData
#' @aliases PrepareData
setMethod("PrepareData", signature(object = "SingleCellExperiment"),
          PrepareData.SingleCellExperiment)


#' @title Principal Component Analysis
#'
#' @description Perform principal component analysis using assays or the joint 
#' probability matrix as input.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param assay.name Name of the assay to compute PCA. One of \code{assayNames(object)}
#' or \code{joint.probability}. By default \code{joint.probability} is used. Use 
#' \code{joint.probability} to obtain an integrated embedding after running 
#' \code{RunParallelDivisiveICP}. One of the assays in \code{assayNames(object)}
#' can be provided before performing integration to assess if data requires 
#' integration.  
#' @param p A positive integer denoting the number of principal components to 
#' calculate and select. Default is \code{50}.
#' @param scale A logical specifying whether the probabilities should be
#' standardized to unit-variance before running PCA. Default is \code{TRUE}.
#' @param center A logical specifying whether the probabilities should be
#' centered before running PCA. Default is \code{TRUE}.
#' @param threshold A threshold for filtering out ICP runs before PCA with the 
#' lower terminal projection accuracy below the threshold. Default is \code{0}.
#' @param pca.method A character specifying the PCA method. One of \code{"irlba"}
#' (default), \code{"RSpectra"} or \code{"stats"}. Set seed before, if the method 
#' is \code{"irlba"} to ensure reproducibility. 
#' @param return.model A logical specifying if the PCA model should or not be 
#' retrieved. By default \code{FALSE}. Only implemented for \code{pca.method = "stats"}. 
#' If \code{TRUE}, the \code{pca.method} is coerced to \code{"stats"}. 
#' @param select.icp.tables Select the ICP cluster probability tables to perform 
#' PCA. By default \code{NULL}, i.e., all are used, except if the ICP tables were
#' obtained with the function \code{RunParallelDivisiveICP}, in which the ICP 
#' tables correspond to the last round of divisive clustering for every epoch.  
#' A vector of \code{integers} should be given otherwise.   
#' @param features A character of feature names matching \code{row.names(object)} 
#' to select from before computing PCA. Only used if \code{assay.name} is one of
#' the assays in \code{assayNames(object)}, otherwise it is ignored. 
#' @param dimred.name Dimensional reduction name given to the returned PCA. By 
#' default \code{"PCA"}. 
#'
#' @name RunPCA
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords PCA eigendecomposition
#'
#' @importFrom RSpectra eigs_sym
#' @importFrom irlba prcomp_irlba
#' @importFrom stats prcomp
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay assayNames
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Prepare SCE object for analysis
#' sce <- PrepareData(sce)
#' 
#' # Multi-level integration (just for highlighting purposes; use default parameters)
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "Batch", 
#'                               k = 2, L = 25, C = 1, train.k.nn = 10, 
#'                               train.k.nn.prop = NULL, use.cluster.seed = FALSE,
#'                               build.train.set = FALSE, ari.cutoff = 0.1, 
#'                              threads = 2)
#' 
#' # Integrated PCA
#' set.seed(125) # to ensure reproducibility for the default 'irlba' method
#' sce <- RunPCA(object = sce, assay.name = "joint.probability", p = 10)
#' 
#' # Plot result 
#' cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
#'                               legend.nrow = 1),
#'                    PlotDimRed(object = sce, color.by = "Species", 
#'                              legend.nrow = 1), ncol = 2)
#'
RunPCA.SingleCellExperiment <- function(object, assay.name, p, scale, center, threshold,
                                        pca.method, return.model, select.icp.tables, 
                                        features, dimred.name) {
    
    # Check arguments
    stopifnot(is(object, "SingleCellExperiment"), 
             (is.character(assay.name) && (length(assay.name)==1) && ( ((assay.name == "joint.probability") && ("joint.probability" %in% names(metadata(object)$coralysis))) || (assay.name %in% assayNames(object)))),
             (is.numeric(p) && (length(p)==1) && (p%%1==0)), 
             is.logical(scale), is.logical(center), 
             ((threshold >= 0) && (threshold < 1)),
             (is.character(pca.method) && (length(pca.method)==1) && (pca.method %in% c("irlba", "RSpectra", "stats"))), 
             is.logical(return.model), 
             (is.null(select.icp.tables) || (is.numeric(select.icp.tables) && all(select.icp.tables %% 1 == 0))), 
             (is.null(features) || (is.character(features) && all(features %in% row.names(object)))), 
             (is.character(dimred.name) && (length(dimred.name)==1)))
    
    # Retrieve cell matrix
    if (assay.name == "joint.probability") { # select probability from 'metadata(object)$coralysis$joint.probability'
        # Get ICP tables
        n.icps <- length(metadata(object)$coralysis$joint.probability)
        if (is.null(select.icp.tables)) {
            select.icp.tables <- 1:n.icps
            divisive.icp <- metadata(object)$coralysis$divisive.icp
            if (!is.null(divisive.icp)) {
                L <- metadata(object)$coralysis$L
                k <- metadata(object)$coralysis$k  
                Ks <- log2(k) # divisive K rounds (if k=16: 2 --> 4 --> 8 --> 16, i.e., Ks = 4 rounds)
                select.icp.tables <- seq(Ks, Ks * L, Ks) # select a ICP table every Ks round
                message(paste0("Divisive ICP: selecting ICP tables multiple of ", Ks))
            }
        }
        if (threshold == 0) {
            X <- do.call(cbind, metadata(object)$coralysis$joint.probability[select.icp.tables])
        } else {
            icp_runs_logical <- unlist(lapply(metadata(object)$coralysis$metrics, function(x) x["ARI",])) >= threshold
            icp_runs_logical <- (icp_runs_logical & ((1:n.icps) %in% select.icp.tables))
            X <- do.call(cbind, metadata(object)$coralysis$joint.probability[icp_runs_logical])
        }
    } else { # select assay from 'assayNames(object)'
        if (is.null(features)) {
            features <- row.names(object)
        } 
        X <- t(assay(object[features,], assay.name))
    }
    
    # Check if assumptions are met
    if (p > ncol(X)) {
        stop(paste0("p larger than number of features. Decrease p"))
    }
    if (return.model) {
        if (pca.method != "stats") {
            message(paste0("Setting 'pca.method' to 'stats' as 'return.model' is TRUE.", 
                           "\nSet 'return.model' to FALSE to use 'method' '", pca.method, "'."))
            pca.method <- "stats"
        }
    }

    # Calculate PCA
    if (pca.method=="RSpectra") {
        X <- scale(X, scale = scale, center = center)
        # X^T %*% X
        A = crossprod(X)
        # Perform eigen decomposition
        eigs_sym_out <- RSpectra::eigs_sym(A, p, which = "LM")
        pca <- X %*% eigs_sym_out$vectors
        colnames(pca) <- paste0("PC", seq_len(ncol(pca))) 
    } 
    if (pca.method == "irlba") {
        pca <- irlba::prcomp_irlba(x = X, n = p, scale. = scale, center = center)
        pca <- pca$x
    }
    if (pca.method == "stats") {
        pca <- stats::prcomp(x = X, scale. = scale, center = center, rank = p)
        if (return.model) {
            metadata(object)$coralysis$pca.model <- pca
        }
        pca <- pca$x
    }
    
    # Saving PCA into SCE object (& params)
    reducedDim(object, type = dimred.name) <- pca
    metadata(object)$coralysis$p <- p # saving due to compatibility issues
    metadata(object)$coralysis$pca.params <-  list("p" = p, 
                                                   "assay.name" = assay.name,
                                                   "scale" = scale, 
                                                   "threshold" = threshold, 
                                                   "center" = center,
                                                   "pca.method" = pca.method, 
                                                   "return.model" = return.model, 
                                                   "select.icp.tables" = select.icp.tables, 
                                                   "dimred.name" = dimred.name)
    
    return(object)
}
#' @rdname RunPCA
#' @aliases RunPCA
setMethod("RunPCA", signature(object = "SingleCellExperiment"),
          RunPCA.SingleCellExperiment)


#' @title Uniform Manifold Approximation and Projection (UMAP)
#'
#' @description Run nonlinear dimensionality reduction using UMAP with a dimensional 
#' reduction as input.
#'
#' @param object An object of \code{SingleCellExperiment} class.
#' @param dims Dimensions to select from \code{dimred.type}. By default \code{NULL}, 
#' i.e., all the dimensions are selected. Provide a numeric vector to select a 
#' specific range, e.g., \code{dims = 1:10} to select the first 10 dimensions. 
#' @param dimred.type Dimensional reduction type to use. By default \code{"PCA"}.
#' @param return.model Return UMAP model. By default \code{FALSE}.
#' @param umap.method UMAP method to use: \code{"umap"} or \code{"uwot"}. 
#' By default \code{"umap"}.
#' @param dimred.name Dimensional reduction name given to the returned UMAP. 
#' By default \code{"UMAP"}.  
#' @param ... Parameters to be passed to the \code{umap} function. The parameters 
#' given should match the parameters accepted by the \code{umap} function depending
#' on the \code{umap.method} given. Check possible parameters with \code{?umap::umap} 
#' or \code{?uwot::umap} depending if \code{umap.method} is \code{"umap"} or 
#' \code{"uwot"}.  
#'
#' @name RunUMAP
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @keywords Uniform Manifold Approximation and Projection UMAP
#'
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Prepare SCE object for analysis
#' sce <- PrepareData(sce)
#' 
#' # Multi-level integration (just for highlighting purposes; use default parameters)
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "Batch", 
#'                               k = 2, L = 25, C = 1, train.k.nn = 10, 
#'                               train.k.nn.prop = NULL, use.cluster.seed = FALSE,
#'                               build.train.set = FALSE, ari.cutoff = 0.1, 
#'                              threads = 2)
#' 
#' # Integrated PCA
#' set.seed(125) # to ensure reproducibility for the default 'irlba' method
#' sce <- RunPCA(object = sce, assay.name = "joint.probability", p = 10)
#' 
#' # Plot result 
#' cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
#'                               legend.nrow = 1),
#'                    PlotDimRed(object = sce, color.by = "Species", 
#'                              legend.nrow = 1), ncol = 2)
#'
#' # Run UMAP
#' set.seed(123)
#' sce <- RunUMAP(sce, dimred.type = "PCA")
#' 
#' # Plot results
#' # Plot result 
#' cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
#'                               legend.nrow = 1),
#'                    PlotDimRed(object = sce, color.by = "Species", 
#'                              legend.nrow = 1), ncol = 2)
#'                              
RunUMAP.SingleCellExperiment <- function(object, dims, dimred.type, return.model, 
                                         umap.method, dimred.name, ...) {
    
    ## Check arguments
    # Check 'return.model', 'dimred.name', 'umap.method' & 'dimred.type'
    stopifnot(is.logical(return.model), is.character(dimred.name), is.character(umap.method), is.character(dimred.type))
    if ( ! (dimred.type %in% reducedDimNames(object)) ) {
        stop(paste("The 'dimred.type' '", dimred.type, "' is not among 'reducedDimNames(object)':", 
                   paste(reducedDimNames(x = object), collapse = ", ")))
    }
    dimred <- reducedDim(x = object, type = dimred.type)
    if (is.null(dims)) { # select all the dimensions from 'dimred.type' if 'dims = NULL'
        dims <- seq_len(ncol(dimred))
    }
    # Check 'umap.method'
    if ( ! (umap.method %in% c("umap", "uwot")) ) {
        stop("'umap.method' must be one of: 'umap', 'uwot'")
    }
    
    # Select UMAP method function
    if (umap.method == "umap") {
        umap.fun <- umap::umap
        slot.name <- "layout"
        params <- c(list(d = dimred[,dims]), ...)
    }
    if (umap.method == "uwot") {
        umap.fun <- uwot::umap
        slot.name <- "embedding"
        params <- c(list(X = dimred[,dims], ret_model = return.model), ...)
    }
    
    # Run UMAP
    umap.out <- do.call(umap.fun, params) 
    
    # Model
    if (return.model) {
        metadata(object)$coralysis$umap.model <- umap.out
    } else {
        if (umap.method == "uwot") { # Turn matrix into a list to match 'umap' output structure 
            umap.out <- list("embedding" = umap.out)
        }
    }
    
    # Return SCE object with UMAP
    reducedDim(x = object, type = dimred.name) <- umap.out[[slot.name]]
    return(object)
}
#' @rdname RunUMAP
#' @aliases RunUMAP
setMethod("RunUMAP", signature(object = "SingleCellExperiment"),
          RunUMAP.SingleCellExperiment)


#' @title Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding 
#' (t-SNE)
#'
#' @description Run nonlinear dimensionality reduction using t-SNE with the 
#' PCA-transformed consensus matrix as input.
#'
#' @param object Object of \code{SingleCellExperiment} class.
#' @param dims Dimensions to select from \code{dimred.type}. By default \code{NULL}, 
#' i.e., all the dimensions are selected. Provide a numeric vector to select a 
#' specific range, e.g., \code{dims = 1:10} to select the first 10 dimensions. 
#' @param dimred.type Dimensional reduction type to use. By default \code{"PCA"}.
#' @param perplexity Perplexity of t-SNE.
#' @param dimred.name Dimensional reduction name given to the returned t-SNE. 
#' By default \code{"TSNE"}.  
#' @param ... Parameters to be passed to the \code{Rtsne} function. The parameters 
#' given should match the parameters accepted by the \code{Rtsne} function. Check 
#' possible parameters with \code{?Rtsne::Rtsne}.
#'
#' @name RunTSNE
#'
#' @return A \code{SingleCellExperiment} object. 
#'
#' @keywords Barnes-Hut implementation of t-Distributed
#' Stochastic Neighbor Embedding t-SNE
#'
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Run PCA
#' set.seed(125) # to ensure reproducibility for the default 'irlba' method
#' sce <- RunPCA(object = sce, assay.name = "logcounts", 
#'               pca.method = "stats", p = nrow(sce))
#' 
#' # Run t-SNE
#' set.seed(125) # to ensure reproducibility for the default 'irlba' method
#' sce <- RunTSNE(object = sce, dimred.type = "PCA", check_duplicates = FALSE)
#' 
#' # Plot result 
#' cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
#'                               legend.nrow = 1),
#'                    PlotDimRed(object = sce, color.by = "Species", 
#'                              legend.nrow = 1), ncol = 2)
#'                    
RunTSNE.SingleCellExperiment <- function(object, dims, dimred.type, perplexity, dimred.name, ...) {
    
    # Check arguments 
    stopifnot(is(object, "SingleCellExperiment"), any(is.null(dims), is.numeric(dims)), is.numeric(perplexity),
              all(is.character(dimred.type), (dimred.type %in% reducedDimNames(object))))
    if ( ! (dimred.type %in% reducedDimNames(object)) ) {
        stop(paste("The 'dimred.type' '", dimred.type, "' is not among 'reducedDimNames(object)':", 
                   paste(reducedDimNames(x = object), collapse = ", ")))
    }
    
    # Parameters
    dimred <- reducedDim(x = object, type = dimred.type)
    if (is.null(dims)) { # select all the dimensions from 'dimred.type' if 'dims = NULL'
        dims <- seq_len(ncol(dimred))
    }
    params <- c(list(X = dimred[,dims], is_distance = FALSE, perplexity = perplexity, pca = FALSE), ...)
    
    # Run t-SNE
    rtsne.out <- do.call(Rtsne::Rtsne, params)
    
    # Return t-SNE
    reducedDim(x = object, type = dimred.name) <- rtsne.out$Y
    return(object)
}
#' @rdname RunTSNE
#' @aliases RunTSNE
setMethod("RunTSNE", signature(object = "SingleCellExperiment"),
          RunTSNE.SingleCellExperiment)


#' @title Identification of feature markers for all clusters
#'
#' @description \code{FindAllClusterMarkers} enables identifying feature markers for all 
#' clusters at once. This is done by differential expresission analysis where 
#' cells from one cluster are compared against the cells from the rest of the 
#' clusters. Feature and cell filters can be applied to accelerate the analysis, 
#' but this might lead to missing weak signals.
#'
#' @param object A \code{SingleCellExperiment} object. 
#' @param clustering.label A variable name (of class \code{character}) available 
#' in the cell metadata \code{colData(object)} with the clustering labels 
#' (\code{character} or \code{factor}) to use.
#' @param test Which test to use. Only "wilcox" (the Wilcoxon rank-sum test,
#' AKA Mann-Whitney U test) is supported at the moment.
#' @param log2fc.threshold Filters out features that have log2 fold-change of the
#' averaged feature expression values below this threshold. Default is \code{0.25}.
#' @param min.pct Filters out features that have dropout rate (fraction of cells
#' expressing a feature) below this threshold in both comparison groups. Default is 
#' \code{0.1}.
#' @param min.diff.pct Filters out features that do not have this minimum
#' difference in the dropout rates (fraction of cells expressing a feature)
#' between the two comparison groups. Default is \code{NULL}.
#' @param min.cells.group The minimum number of cells in the two comparison
#' groups to perform the DE analysis. If the number of cells is below the
#' threshold, then the DE analysis of this cluster is skipped. Default is \code{3}.
#' @param max.cells.per.cluster The maximum number of cells per cluster if
#' downsampling is performed to speed up the DE analysis. Default is \code{NULL}, 
#' i.e., no downsampling.
#' @param return.thresh If \code{only.pos=TRUE}, then return only features that have the
#' adjusted p-value (adjusted by the Bonferroni method) below or equal to this
#' threshold. Default is \code{0.01}.
#' @param only.pos Whether to return only features that have an adjusted p-value
#' (adjusted by the Bonferroni method) below or equal to the threshold. Default 
#' is \code{FALSE}.
#'
#' @name FindAllClusterMarkers
#'
#' @return A data frame of the results if positive results were found, else \code{NULL}.
#'
#' @keywords differential expression DE analysis feature markers
#'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment logcounts
#' @importFrom stats wilcox.test p.adjust
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Markers
#' dge <- FindAllClusterMarkers(sce, clustering.label = "Species")
#' dge
#' 
FindAllClusterMarkers.SingleCellExperiment <- function(object,
                                                       clustering.label,
                                                       test,
                                                       log2fc.threshold,
                                                       min.pct,
                                                       min.diff.pct,
                                                       min.cells.group,
                                                       max.cells.per.cluster,
                                                       return.thresh,
                                                       only.pos) {
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), 
              all(is.character(clustering.label), length(clustering.label)==1, clustering.label %in% colnames(colData(object)), (is.character(object[[clustering.label]]) || is.factor(object[[clustering.label]]))), 
              test == "wilcox", 
              (is.numeric(log2fc.threshold) && (length(log2fc.threshold)==1) && (log2fc.threshold > 0)), 
              (is.numeric(min.pct) && (length(min.pct)==1) && ((min.pct >= 0) && (min.pct <= 1))), 
              (is.null(min.diff.pct) || (is.numeric(min.diff.pct) && (length(min.diff.pct)==1) && ((min.diff.pct >= 0) && (min.diff.pct <= 1)))), 
              (is.numeric(min.cells.group) && (length(min.cells.group)==1) && (min.cells.group > 0)), 
              (is.null(max.cells.per.cluster) || (is.numeric(max.cells.per.cluster) && (length(max.cells.per.cluster)==1) && (max.cells.per.cluster > 0))), 
              (is.numeric(return.thresh) && (length(return.thresh)==1) && (return.thresh > 0)), 
              is.logical(only.pos) && (length(only.pos)==1))
    
    # Retrieve data & no. of features
    number.of.expressed.features <- nrow(object)
    data <- logcounts(object)
    
    # Retrieve clustering as factor
    clustering <- object[[clustering.label]]
    if (is.character(clustering)) { # if character, coerce to factor
        clustering <- as.factor(clustering)
    }
    names(clustering) <- colnames(object)
    clusters <- levels(clustering)
    
    # Downsample each cluster
    if (!is.null(max.cells.per.cluster)) {
        cells_downsampled <- c()
        for (cluster in clusters) {
            cells_in_cluster <- table(clustering)[cluster]
            if (max.cells.per.cluster < cells_in_cluster) {
                inds <- sample(seq_len(cells_in_cluster), size = max.cells.per.cluster, replace = FALSE)
                names_cluster <- names(clustering[clustering==cluster])
                cells_downsampled <- c(cells_downsampled,names_cluster[inds])
            }
        }
        data <- data[,cells_downsampled]
        clustering <- clustering[cells_downsampled]
    }
    
    # Compare cells from each cluster against all other clusters
    results_list <- list()
    for (cluster in clusters) {
        cat("-----------------------------------\n")
        cat(paste0("testing cluster ",cluster,"\n"))
        
        # Extract data
        data_cluster <- data[,clustering==cluster]
        data_other <- data[,clustering!=cluster]
        
        # Skip if the number of cells in the test
        # or the reference set is lower than min.cells.group
        if (ncol(data_cluster) < min.cells.group | ncol(data_other) < min.cells.group) {
            cat("-----------------------------------\n")
            next
        }
        # min.pct filter
        features.pct_cluster <- apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
        features.pct_other <- apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)
        features_to_include <- rownames(data_cluster)[features.pct_cluster>=min.pct | features.pct_other >= min.pct]
        data_cluster <- data_cluster[features_to_include,,drop=FALSE]
        data_other <- data_other[features_to_include,,drop=FALSE]
        cat(paste0(nrow(data_cluster)," features left after min.pct filtering\n"))
        if (nrow(data_cluster)==0) {
            cat("-----------------------------------\n")
            next
        }
        
        # min.diff.pct filter
        if (!is.null(min.diff.pct)) {
            features.pct_cluster <- features.pct_cluster[features_to_include]
            features.pct_other <- features.pct_other[features_to_include]
            features_to_include <- rownames(data_cluster)[abs(features.pct_cluster-features.pct_other) >= min.diff.pct]
            data_cluster <- data_cluster[features_to_include,,drop=FALSE]
            data_other <- data_other[features_to_include,,drop=FALSE]
        }
        cat(paste0(nrow(data_cluster)," features left after min.diff.pct filtering\n"))
        if (nrow(data_cluster)==0) {
            cat("-----------------------------------\n")
            next
        }
        
        # logfc.threshold filter
        # Calculate log2 fold changes
        cluster_aves <- apply(data_cluster,1,mean)
        other_aves <- apply(data_other,1,mean)
        log2FC <- cluster_aves - other_aves
        features_to_include <- rownames(data_cluster)[log2FC >= log2fc.threshold | log2FC <= -log2fc.threshold]
        data_cluster <- data_cluster[features_to_include,,drop=FALSE]
        data_other <- data_other[features_to_include,,drop=FALSE]
        cat(paste0(nrow(data_cluster)," features left after log2fc.threshold filtering\n"))
        if (nrow(data_cluster)==0) {
            cat("-----------------------------------\n")
            next
        }
        
        # Run DE test
        if (test=="wilcox") {
            wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,], y=data_other[x,]))
            p_values <- unlist(lapply(wilcox.res, function(x) x$p.value))
            names(p_values) <- rownames(data_cluster)
            
            # Adjust p-values
            adj_p_values <- p.adjust(p_values, method = "bonferroni", n = number.of.expressed.features)
            
            res <- cbind(p_values,adj_p_values,log2FC[names(p_values)],features.pct_cluster[names(p_values)],features.pct_other[names(p_values)],abs(features.pct_cluster-features.pct_other)[names(p_values)])
            colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
            res <- as.data.frame(res)
            res$cluster <- cluster
            res$marker <- names(p_values)
        }
        results_list[[cluster]] <- res
        cat("-----------------------------------\n")
    }
    
    results_df <- do.call(rbind,results_list)
    rownames(results_df) <- make.unique(unlist(lapply(results_list,rownames)))
    
    if(only.pos) {
        results_df <- results_df[results_df$adj.p.value <= return.thresh,]
        return(results_df)
    }
    
    return(results_df)
}
#' @rdname FindAllClusterMarkers
#' @aliases FindAllClusterMarkers
setMethod("FindAllClusterMarkers", signature(object = "SingleCellExperiment"),
          FindAllClusterMarkers.SingleCellExperiment)


#' @title Differential expression between cell clusters
#'
#' @description \code{FindClusterMarkers} enables identifying feature markers for one 
#' cluster or two arbitrary combinations of clusters, e.g. 1_2 vs. 3_4_5. Feature 
#' and cell filters can be applied to accelerate the analysis, but this might 
#' lead to missing weak signals.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param clustering.label A variable name (of class \code{character}) available 
#' in the cell metadata \code{colData(object)} with the clustering labels 
#' (\code{character} or \code{factor}) to use.
#' @param clusters.1 a character or numeric vector denoting which clusters
#' to use in the first group (named group.1 in the results)
#' @param clusters.2 a character or numeric vector denoting which clusters
#' to use in the second group (named group.2 in the results)
#' @param test Which test to use. Only "wilcoxon" (the Wilcoxon rank-sum test,
#' AKA Mann-Whitney U test) is supported at the moment.
#' @param log2fc.threshold Filters out features that have log2 fold-change of the
#' averaged feature expression values below this threshold.
#' Default is \code{0.25}.
#' @param min.pct Filters out features that have dropout rate (fraction of cells
#' expressing a feature) below this threshold in both comparison groups
#' Default is \code{0.1}.
#' @param min.diff.pct Filters out features that do not have this minimum
#' difference in the dropout rates (fraction of cells expressing a feature)
#' between the two comparison groups. Default is \code{NULL}.
#' @param min.cells.group The minimum number of cells in the two comparison
#' groups to perform the DE analysis. If the number of cells is below the
#' threshold, then the DE analysis is not performed.
#' Default is \code{3}.
#' @param max.cells.per.cluster The maximun number of cells per cluster
#' if downsampling is performed to speed up the DE analysis.
#' Default is \code{NULL}, i.e. no downsampling.
#' @param return.thresh If only.pos=TRUE, then return only features that
#' have the adjusted p-value (adjusted by the Bonferroni method) below or
#' equal to this threshold.  Default is \code{0.01}.
#' @param only.pos Whether to return only features that have an adjusted
#' p-value (adjusted by the Bonferroni method) below or equal to the
#' threshold. Default is \code{FALSE}.
#'
#' @name FindClusterMarkers
#'
#' @return a data frame of the results if positive results were found, else NULL
#'
#' @keywords differential expression DE analysis feature markers
#'
#' @importFrom S4Vectors metadata
#' @importFrom stats wilcox.test p.adjust
#' @importFrom SingleCellExperiment logcounts
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Markers between versicolor vs virginica
#' dge <- FindClusterMarkers(sce, clustering.label = "Species", 
#'                           clusters.1 = "versicolor", 
#'                           clusters.2 = "virginica")
#' dge
#' 
FindClusterMarkers.SingleCellExperiment <- function(object, 
                                                    clustering.label,
                                                    clusters.1,
                                                    clusters.2,
                                                    test,
                                                    log2fc.threshold,
                                                    min.pct,
                                                    min.diff.pct,
                                                    min.cells.group,
                                                    max.cells.per.cluster,
                                                    return.thresh,
                                                    only.pos) {
    
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), 
              all(is.character(clustering.label), length(clustering.label)==1, clustering.label %in% colnames(colData(object)), (is.character(object[[clustering.label]]) || is.factor(object[[clustering.label]]))), 
              all(
                  (is.null(clusters.1) || (is.character(clusters.1) && (length(clusters.1)==1) && (clusters.1 %in% object[[clustering.label]]))) && 
                      (is.null(clusters.2) || (is.character(clusters.2) && (length(clusters.2)==1) && (clusters.2 %in% object[[clustering.label]])))
              ),
              test == "wilcox", 
              (is.numeric(log2fc.threshold) && (length(log2fc.threshold)==1) && (log2fc.threshold > 0)), 
              (is.numeric(min.pct) && (length(min.pct)==1) && ((min.pct >= 0) && (min.pct <= 1))), 
              (is.null(min.diff.pct) || (is.numeric(min.diff.pct) && (length(min.diff.pct)==1) && ((min.diff.pct >= 0) && (min.diff.pct <= 1)))), 
              (is.numeric(min.cells.group) && (length(min.cells.group)==1) && (min.cells.group > 0)), 
              (is.null(max.cells.per.cluster) || (is.numeric(max.cells.per.cluster) && (length(max.cells.per.cluster)==1) && (max.cells.per.cluster > 0))), 
              (is.numeric(return.thresh) && (length(return.thresh)==1) && (return.thresh > 0)), 
              is.logical(only.pos) && (length(only.pos)==1))
    
    # Retrieve clustering as factor
    clustering <- object[[clustering.label]]
    if (is.character(clustering)) { # if character, coerce to factor
        clustering <- as.factor(clustering)
    }
    names(clustering) <- colnames(object)
    clusters <- levels(clustering)
    
    # Retrieve data & cells
    data <- logcounts(object)
    cells_to_include_1 <- names(clustering)[clustering %in% clusters.1]
    clustering_1 <- factor(rep("group.1",length(cells_to_include_1)))
    names(clustering_1) <- cells_to_include_1
    if (is.null(clusters.2)) {
        clusters.2 <- setdiff(levels(clustering),clusters.1)
    }
    cells_to_include_2 <- names(clustering)[clustering %in% clusters.2]
    clustering_2 <- factor(rep("group.2",length(cells_to_include_2)))
    names(clustering_2) <- cells_to_include_2
    data <- data[,c(cells_to_include_1,cells_to_include_2)]
    clustering <- factor(c(as.character(clustering_1),as.character(clustering_2)))
    names(clustering) <- c(cells_to_include_1,cells_to_include_2)
    clusters <- levels(clustering)
    
    # Remove features that are not expressed in any of the cells
    data <- data[Matrix::rowSums(data)!=0,]
    clusters <- levels(clustering)
    
    # Downsample each cluster
    if (!is.null(max.cells.per.cluster)) {
        cells_downsampled <- c()
        for (cluster in clusters) {
            cells_in_cluster <- table(clustering)[cluster]
            if (max.cells.per.cluster < cells_in_cluster) {
                inds <- sample(seq_len(cells_in_cluster),
                               size = max.cells.per.cluster,
                               replace = FALSE)
                cells_downsampled <- c(cells_downsampled,
                                       names(clustering[clustering==cluster])[inds])
            }
        }
        data <- data[,cells_downsampled]
        clustering <- clustering[cells_downsampled]
    }
    
    # Compare cells from each cluster against all other clusters
    results_list <- list()
    cluster <- "group.1"
    cat(paste0("testing cluster ",cluster,"\n"))
    
    # Extract data
    data_cluster <- data[,clustering==cluster]
    data_other <- data[,clustering!=cluster]
    
    # Skip if the number of cells in the test or the reference set
    # is lower than min.cells.group
    if (ncol(data_cluster) < min.cells.group | ncol(data_other) < min.cells.group) {
        cat("-----------------------------------\n")
        return(NULL)
    }
    
    # min.pct filter
    features.pct_cluster <- apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
    features.pct_other <- apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)
    features_to_include <- rownames(data_cluster)[features.pct_cluster>=min.pct | features.pct_other >= min.pct]
    data_cluster <- data_cluster[features_to_include,,drop=FALSE]
    data_other <- data_other[features_to_include,,drop=FALSE]
    cat(paste0(nrow(data_cluster)," features left after min.pct filtering\n"))
    if (nrow(data_cluster)==0) {
        cat("-----------------------------------\n")
        return(NULL)
    }
    
    # min.diff.pct filter
    if (!is.null(min.diff.pct)) {
        features.pct_cluster <- features.pct_cluster[features_to_include]
        features.pct_other <- features.pct_other[features_to_include]
        features_to_include <- rownames(data_cluster)[abs(features.pct_cluster-features.pct_other) >= min.diff.pct]
        data_cluster <- data_cluster[features_to_include,,drop=FALSE]
        data_other <- data_other[features_to_include,,drop=FALSE]
    }
    
    cat(paste0(nrow(data_cluster)," features left after min.diff.pct filtering\n"))
    if (nrow(data_cluster)==0) {
        cat("-----------------------------------\n")
        return(NULL)
    }
    
    # log2fc.threshold filter
    # Calculate log2 fold changes
    cluster_aves <- apply(data_cluster,1,mean)
    other_aves <- apply(data_other,1,mean)
    log2FC <- cluster_aves - other_aves
    features_to_include <- rownames(data_cluster)[log2FC >= log2fc.threshold | log2FC <= -log2fc.threshold]
    data_cluster <- data_cluster[features_to_include,,drop=FALSE]
    data_other <- data_other[features_to_include,,drop=FALSE]
    cat(paste0(nrow(data_cluster)," features left after log2fc.threshold filtering\n"))
    if (nrow(data_cluster)==0) {
        cat("-----------------------------------\n")
        return(NULL)
    }
    
    # Run DE test
    if (test=="wilcox") {
        wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,],y=data_other[x,]))
        p_values <- unlist(lapply(wilcox.res,function(x) x$p.value))
        names(p_values) <- rownames(data_cluster)
        
        # Adjust p-values
        adj_p_values <- p.adjust(p_values, method = "bonferroni", n = nrow(object))
        res <- cbind(p_values,
                     adj_p_values,
                     log2FC[names(p_values)],
                     features.pct_cluster[names(p_values)],
                     features.pct_other[names(p_values)],
                     abs(features.pct_cluster-features.pct_other)[names(p_values)])
        colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
        res <- as.data.frame(res)
        res$cluster <- cluster
        res$marker <- names(p_values)
    }
    
    results_list[[cluster]] <- res
    results_df <- do.call(rbind,results_list)
    rownames(results_df) <- make.unique(unlist(lapply(results_list,rownames)))
    results_df$cluster <- NULL
    if(only.pos) {
        results_df <- results_df[results_df$adj.p.value <= return.thresh,]
        return(results_df)
    }
    
    return(results_df)
}
#' @rdname FindClusterMarkers
#' @aliases FindClusterMarkers
setMethod("FindClusterMarkers", signature(object = "SingleCellExperiment"),
          FindClusterMarkers.SingleCellExperiment)


#' @title Get ICP cell cluster probability
#' 
#' @description Get ICP cell cluster probability table(s)
#' 
#' @param object An object of \code{SingleCellExperiment} class with ICP cell 
#' cluster probability tables saved in \code{metadata(object)$coralysis$joint.probability}. 
#' After running \code{RunParallelDivisiveICP}. 
#' @param icp.run ICP run(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Specify a numeric vector to 
#' retrieve a specific set of tables. 
#' @param icp.round ICP round(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. 
#' @param concatenate Concatenate list of ICP cell cluster probability tables retrieved. 
#' By default \code{TRUE}, i.e., the list of ICP cell cluster probability tables is
#' concatenated. 
#' 
#' @name GetCellClusterProbability
#' 
#' @return A list with ICP cell cluster probability tables or a matrix with 
#' concatenated tables.  
#' 
#' @keywords Cell cluster probability
#'  
#' @importFrom S4Vectors metadata metadata<-
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Prepare SCE object for analysis
#' sce <- PrepareData(sce)
#' 
#' # Multi-level integration (just for highlighting purposes; use default parameters)
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "Batch", 
#'                               k = 2, L = 25, C = 1, train.k.nn = 10, 
#'                               train.k.nn.prop = NULL, use.cluster.seed = FALSE,
#'                               build.train.set = FALSE, ari.cutoff = 0.1, 
#'                              threads = 2)
#' 
#' # Get cluster probability for all ICP runs 
#' probs <- GetCellClusterProbability(object = sce, icp.round = 1, concatenate = TRUE) 
#' probs[1:10, 1:5]
#' 
GetCellClusterProbability.SingleCellExperiment <- function(object, icp.run, icp.round, concatenate) {
    
    # Retrieve important params
    L <- metadata(object)$coralysis$L
    k <- metadata(object)$coralysis$k
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), is(metadata(object)$coralysis$joint.probability, "list"), 
              is.numeric(L), is.numeric(k), any(is.null(icp.run), (is.numeric(icp.run) && all(icp.run <= L))), 
              any(is.null(icp.round), (is.numeric(icp.round))), is.logical(concatenate))
    
    # If ICP run is NULL retrieve all
    if (is.null(icp.run)) {
        icp.run <- seq_len(L)
    }
    
    # If divisive ICP, retrieve the right ICP run round
    divisive.icp <- metadata(object)$coralysis$divisive.icp
    if (isTRUE(divisive.icp)) { 
        rounds <- log2(k)
        stopifnot(all(icp.round <= rounds)) # all icp.round need to be equal or lower than rounds
        if (is.null(icp.round)) {
            icp.round <- seq_len(rounds)
        }
        icp.run.tbl <- matrix(data = seq_len(L*rounds), ncol = rounds, byrow = TRUE)
        pick.icp <- c(t(icp.run.tbl[icp.run, icp.round]))
    } else { # if not divisive just select icp runs
        pick.icp <- icp.run
    }
    
    # Return probabilities
    probs <- metadata(object)$coralysis$joint.probability[pick.icp]
    if (concatenate) {
        probs <- do.call(cbind, probs)
    }
    return(probs)
}
#' @rdname GetCellClusterProbability
#' @aliases GetCellClusterProbability
setMethod("GetCellClusterProbability", signature(object = "SingleCellExperiment"),
          GetCellClusterProbability.SingleCellExperiment)


#' @title Summarise ICP cell cluster probability
#' 
#' @description Summarise ICP cell cluster probability table(s)
#' 
#' @param object An object of \code{SingleCellExperiment} class with ICP cell 
#' cluster probability tables saved in \code{metadata(object)$coralysis$joint.probability}. 
#' After running \code{RunParallelDivisiveICP}. 
#' @param icp.run ICP run(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Specify a numeric vector to 
#' retrieve a specific set of tables. 
#' @param icp.round ICP round(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. 
#' @param funs Functions to summarise ICP cell cluster probability: \code{"mean"} 
#' and/or \code{"median"}. By default \code{c("mean", "median")}, i.e, both mean
#' and median are calculated. Set to \code{NULL} to not estimate any. 
#' @param scale.funs Scale in the range 0-1 the summarised probability obtained with 
#' \code{funs}. By default \code{TRUE}, i.e., summarised probability will be scaled
#' in the 0-1 range. 
#' @param save.in.sce Save the data frame into the cell metadata from the 
#' \code{SingleCellExperiment} object or return the data frame. By default \code{TRUE}, 
#' i.e., the summary of probabilities retrieved is save in the SCE object in 
#' \code{colData(object)}.  
#' 
#' @name SummariseCellClusterProbability
#' 
#' @return A data frame or a SingleCellExperiment object with ICP cell cluster probability summarised.  
#' 
#' @keywords Summarise cell cluster probability
#' 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
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
#'                                                "Batch" = batch))
#' colnames(sce) <- paste0("samp", 1:ncol(sce))
#' 
#' # Prepare SCE object for analysis
#' sce <- PrepareData(sce)
#' 
#' # Multi-level integration (just for highlighting purposes; use default parameters)
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "Batch", 
#'                               k = 2, L = 25, C = 1, train.k.nn = 10, 
#'                               train.k.nn.prop = NULL, use.cluster.seed = FALSE,
#'                               build.train.set = FALSE, ari.cutoff = 0.1, 
#'                              threads = 2)
#'                              
#' # Integrated PCA
#' set.seed(125) # to ensure reproducibility for the default 'irlba' method
#' sce <- RunPCA(object = sce, assay.name = "joint.probability", p = 10)
#' 
#' # Summarise cluster probability
#' sce <- SummariseCellClusterProbability(object = sce, icp.round = 1, 
#'                                        save.in.sce = TRUE) # saved in 'colData'
#' 
#' # Plot the clustering result for ICP run no. 3 
#' PlotDimRed(object = sce, color.by = "icp_run_round_3_1_clusters")
#' 
#' # Plot Coralysis mean cell cluster probabilities 
#' PlotExpression(object = sce, color.by = "mean_probs", 
#'                color.scale = "viridis")
#' 
SummariseCellClusterProbability.SingleCellExperiment <- function(object, icp.run, icp.round, funs, scale.funs, save.in.sce) {
    
    # Check params that will not be checked in the function below
    stopifnot(any(is.null(funs), any(funs %in% c("mean", "median"))), is.logical(scale.funs), is.logical(save.in.sce))
    
    # Retrieve probability & clustering
    probs <- GetCellClusterProbability(object = object, icp.run = icp.run, icp.round = icp.round, concatenate = FALSE)
    
    # Retrieve important params
    L <- metadata(object)$coralysis$L
    k <- metadata(object)$coralysis$k
    
    # If ICP run is NULL retrieve all
    if (is.null(icp.run)) {
        icp.run <- seq_len(L)
    }
    # If divisive ICP, retrieve the right ICP run round
    divisive.icp <- metadata(object)$coralysis$divisive.icp
    if (isTRUE(divisive.icp)) { 
        rounds <- log2(k)
        if (is.null(icp.round)) {
            icp.round <- seq_len(rounds)
        }
    } else { # if not divisive just select icp runs
        if (is.null(icp.round)) {
            icp.round <- 1
        }
    }
    
    # Summarise
    all.col.names <- c(t(outer(icp.run, icp.round, paste, sep = "_")))
    out <- data.frame(row.names = colnames(object))
    for (p in seq_along(probs)) {
        idx <- apply(X = probs[[p]], MARGIN = 1, FUN = function(x) which.max(x))
        clts <- sapply(X = idx, FUN = function(x) colnames(probs[[p]])[x])
        clt.prob <- apply(X = probs[[p]], MARGIN = 1, FUN = function(x) max(x))
        col.names <- paste("icp_run_round", all.col.names[p], c("clusters", "probs"), sep = "_")
        tmp.tbl <- data.frame(clts, clt.prob)
        colnames(tmp.tbl) <- col.names
        out <- cbind(out, tmp.tbl)
    }
    if ("mean" %in% funs) {
        out$mean_probs <- matrixStats::rowMeans2(as.matrix(out[,grepl("icp_run_round_\\w+_probs", colnames(out))]))
        if (scale.funs) {
            out$scaled_mean_probs <- (((out$mean_probs) - min(out$mean_probs)) / (max(out$mean_probs) - min(out$mean_probs)))
        }
    }
    if ("median" %in% funs) {
        out$median_probs <- matrixStats::rowMedians(as.matrix(out[,grepl("icp_run_round_\\w+_probs", colnames(out))]))
        if (scale.funs) {
            out$scaled_median_probs <- (((out$median_probs) - min(out$median_probs)) / (max(out$median_probs) - min(out$median_probs)))
        }
    }
    
    # Return table
    if (save.in.sce) {
        colData(object) <- cbind(colData(object), out)
        return(object)
    } else{
        return(out)
    }
}
#' @rdname SummariseCellClusterProbability
#' @aliases SummariseCellClusterProbability
setMethod("SummariseCellClusterProbability", signature(object = "SingleCellExperiment"),
          SummariseCellClusterProbability.SingleCellExperiment)


#' @title Get feature coefficients
#' 
#' @description Get feature coefficients from ICP models. 
#' 
#' @param object An object of \code{SingleCellExperiment} class with ICP cell 
#' cluster probability tables saved in \code{metadata(object)$coralysis$joint.probability}. 
#' After running \code{RunParallelDivisiveICP}. 
#' @param icp.run ICP run(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. Specify a numeric vector to 
#' retrieve a specific set of tables. 
#' @param icp.round ICP round(s) to retrieve from \code{metadata(object)$coralysis$joint.probability}. 
#' By default \code{NULL}, i.e., all are retrieved. 
#' 
#' @name GetFeatureCoefficients
#' 
#' @return A list of feature coefficient weights per cluster per ICP run/round.  
#' 
#' @keywords Feature coefficients weights
#'  
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dplyr %>%
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
#' 
#' # GetFeatureCoefficients
#' gene_coefficients_icp_7_1 <- GetFeatureCoefficients(object = sce, icp.run = 7, icp.round = 1)
#' head(gene_coefficients_icp_7_1$icp_13)
#' 
GetFeatureCoefficients.SingleCellExperiment <- function(object, icp.run = NULL, icp.round = NULL) {
    
    # Retrieve important params
    L <- metadata(object)$coralysis$L
    k <- metadata(object)$coralysis$k
    
    # Check input 
    stopifnot(is(object, "SingleCellExperiment"), any(is.null(icp.run), (is.numeric(icp.run) && all(icp.run <= L))), 
              is.numeric(L), is.numeric(k), any(is.null(icp.round), is.numeric(icp.round)))
    
    # Retrieve coefficients
    # If ICP run is NULL retrieve all
    if (is.null(icp.run)) {
        icp.run <- seq_len(L)
    }
    # If divisive ICP, retrieve the right ICP run round
    divisive.icp <- metadata(object)$coralysis$divisive.icp
    if (isTRUE(divisive.icp)) { 
        rounds <- log2(k)
        if (is.null(icp.round)) {
            icp.round <- seq_len(rounds)
        }
        icp.run.tbl <- matrix(data = seq_len(L*rounds), ncol = rounds, byrow = TRUE)
        pick.icp <- c(t(icp.run.tbl[icp.run, icp.round]))
    } else { # if not divisive just select icp runs
        pick.icp <- icp.run
    }
    models <- metadata(object)$coralysis$models[pick.icp]
    feature2coeff <- row.names(object) 
    names(feature2coeff) <- paste0("W", 1:length(feature2coeff))
    feature.coeffs <- lapply(X = models, FUN = function(x) {
        apply(X = x$W[,-ncol(x$W), drop=FALSE], MARGIN = 1, FUN = function(y) {
            idx <- which(y != 0)
            data.frame("feature" = feature2coeff[idx], "coeff" = y[idx])
        })
    }) 
    df.feature.coeffs <- lapply(X = feature.coeffs, FUN = function(x) {
        df <- Reduce(function(y1, y2) merge(y1, y2, by = "feature", all.x = TRUE, all.y = TRUE), x)
        colnames(df) <- c("feature", paste0("coeff_clt", names(x)))
        df[is.na(df)] <- 0
        return(df)
    })
    names(df.feature.coeffs) <- paste0("icp_", pick.icp)
    
    # Return
    return(df.feature.coeffs)
}
#' @rdname GetFeatureCoefficients
#' @aliases GetFeatureCoefficients
setMethod("GetFeatureCoefficients", signature(object = "SingleCellExperiment"),
          GetFeatureCoefficients.SingleCellExperiment)


#' @title Majority voting features by label
#' 
#' @description Get ICP feature coefficients for a label of interest by majority voting label across ICP clusters. 
#' 
#' @param object An object of \code{SingleCellExperiment} class with ICP cell 
#' cluster probability tables saved in \code{metadata(object)$coralysis$joint.probability}. 
#' After running \code{RunParallelDivisiveICP}. 
#' @param label Label of interest available in \code{colData(object)}. 
#' 
#' @name MajorityVotingFeatures
#' 
#' @return A list of with a list of data frames with feature weights per label and a data frame with a summary by label. 
#' 
#' @keywords Majority voting feature coefficients weights
#' 
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom dplyr %>% select filter all_of
#' 
#' @examples 
#' \dontrun{
#' # Import package
#' suppressPackageStartupMessages(library("SingleCellExperiment"))
#' 
#' # Import data from Zenodo
#' data.url <- "https://zenodo.org/records/14845751/files/pbmc_10Xassays.rds?download=1"
#' sce <- readRDS(file = url(data.url))
#' 
#' # Multi-level integration (just for highlighting purposes; use default parameters)
#' set.seed(123)
#' sce <- RunParallelDivisiveICP(object = sce, batch.label = "batch", 
#'                               k = 4, L = 10, C = 1, d = 0.5, 
#'                               train.with.bnn = FALSE, use.cluster.seed = FALSE,
#'                               build.train.set = FALSE, ari.cutoff = 0.1, 
#'                               threads = 2)
#' 
#' # Get coefficients by majority voting for a given categorical variable
#' coeff <- MajorityVotingFeatures(object = sce, label = "cell_type")
#' gene_coeff$summary
#' order.rows <- order(coeff$feature_coeff$Monocyte$coeff_clt2, 
#'                     decreasing = TRUE)
#' head(coeff$feature_coeff$Monocyte[order.rows,], n = 10)
#' }
#' 
MajorityVotingFeatures.SingleCellExperiment <- function(object, label) {
    
    # Check input 
    stopifnot(is(object, "SingleCellExperiment"), (label %in% colnames(colData(object))))
    
    # Cell clusterings
    cell.clts <- SummariseCellClusterProbability(object = object, icp.run = NULL, 
                                                 icp.round = NULL, funs = NULL, 
                                                 scale.funs = FALSE, save.in.sce = FALSE)
    cell.clts <- cell.clts[,grepl("icp_run_round_\\w+_clusters", colnames(cell.clts))]
    
    # Split cell cluster table by label
    cell.clts.label <- split(x = cell.clts, f = colData(object)[, label, drop=TRUE])
    cell.clts.counts <- apply(X = cell.clts, MARGIN = 2, FUN = function(x) table(x))
    clts.label.counts <- lapply(X = cell.clts.label, FUN = function(x) {
        apply(X = x, MARGIN = 2, FUN = function(y) table(y)) 
    })
    
    # geometric mean formula
    geometric_mean <- function(x) exp(mean(log(x)))
    labels <- names(clts.label.counts)
    out <- data.frame("label" = labels, "icp_run" = 0, "icp_round" = 0, "cluster" = "", score = 0)
    res <- clts.label.score <- list()
    for (cell in names(clts.label.counts)) {
        clts.label.score[[cell]] <- list()
        for (icp in names(cell.clts.counts)) {
            tmp <- clts.label.counts[[cell]][[icp]]
            pick.clts <- names(tmp)
            label.per.clts <- (tmp / sum(tmp)) 
            cluster.label <- (tmp / cell.clts.counts[[icp]][pick.clts]) 
            clts.label.score[[cell]][[icp]] <- apply(X = data.frame(label.per.clts, cluster.label)[,c(2,4)], MARGIN = 1, FUN = function(x) geometric_mean(x))
        }
        idx.max <- which.max(unlist(clts.label.score[[cell]]))
        score <- unlist(clts.label.score[[cell]])[idx.max]
        id <- strsplit(x = names(unlist(clts.label.counts[[cell]])[idx.max]), split = "_")[[1]]
        out[out$label==cell, c("icp_run", "icp_round", "score")] <- c(as.numeric(id[4:5]), score)
        out[out$label==cell, "cluster"] <- gsub(pattern = "clusters.", replacement = "", x = id[6])
        # Get coefficients
        pick.clt <- paste0("coeff_clt", out[out$label==cell, "cluster"])
        res[[cell]] <- GetFeatureCoefficients(object = object, icp.run = out[out$label==cell, "icp_run"], 
                                              icp.round = out[out$label==cell, "icp_round"])[[1]] %>% 
            select(all_of(c("feature", pick.clt))) %>% 
            filter(.data[[pick.clt]] != 0)
    }
    
    # Return
    res.out <- list("feature_coeff" = res, "summary" = out)
    return(res.out)
}
#' @rdname MajorityVotingFeatures
#' @aliases MajorityVotingFeatures
setMethod("MajorityVotingFeatures", signature(object = "SingleCellExperiment"),
          MajorityVotingFeatures.SingleCellExperiment)
