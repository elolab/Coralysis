#' @title Multi-level integration
#'
#' @description Run divisive ICP clustering in parallel in order to perform multi-level integration.
#'
#' @param object An object of \code{SingleCellExperiment} class.
#' @param batch.label A variable name (of class \code{character}) available 
#' in the cell metadata \code{colData(object)} with the batch labels (\code{character} 
#' or \code{factor}) to use. The variable provided must not contain \code{NAs}.
#' By default \code{NULL}, i.e., cells are sampled evenly regardless their batch. 
#' @param k A positive integer power of two, i.e., \code{2**n}, where \code{n>0},
#' specifying the number of clusters in the last Iterative Clustering Projection (ICP)
#' round. Decreasing \code{k} leads to smaller cell populations diversity and vice versa. 
#' Default is \code{16}, i.e., the divisive clustering 2 -> 4 -> 8 -> 16 is performed. 
#' @param d A numeric greater than \code{0} and smaller than \code{1} that
#' determines how many cells \code{n} are down- or oversampled from each cluster
#' into the training data (\code{n=N/k*d}), where \code{N} is the total number
#' of cells, \code{k} is the number of clusters in ICP. Increasing above 0.3
#' leads greadually to smaller cell populations diversity.
#' Default is \code{0.3}.
#' @param L A positive integer greater than \code{1} denoting the number of
#' the ICP runs to run. Default is \code{50}. Increasing recommended with
#' a significantly larger sample size (tens of thousands of cells).
#' Default is \code{200}.
#' @param r A positive integer that denotes the number of reiterations
#' performed until the ICP algorithm stops.
#' Increasing recommended with a significantly larger sample size
#' (tens of thousands of cells). Default is \code{5}.
#' @param C A positive real number denoting the cost of constraints violation in
#' the L1-regularized logistic regression model from the LIBLINEAR library.
#' Decreasing leads to more stringent feature selection, i.e. less features are
#' selected that are used to build the projection classifier. Decreasing to a
#' very low value (~ \code{0.01}) can lead to failure to identify central cell
#' populations. Default \code{0.3}.
#' @param reg.type "L1" or "L2". L2-regularization was not
#' investigated in the manuscript, but it leads to a more conventional
#' outcome (less subpopulations). Default is "L1".
#' @param max.iter A positive integer that denotes
#' the maximum number of iterations performed until ICP stops. This parameter
#' is only useful in situations where ICP converges extremely slowly, preventing
#' the algorithm to run too long. In most cases, reaching
#' the number of reiterations (\code{r=5}) terminates the algorithm.
#' Default is \code{200}.
#' @param threads A positive integer that specifies how many logical processors
#' (threads) to use in parallel computation.
#' Set \code{1} to disable parallelism altogether or \code{0} to use all
#' available threas except one. Default is \code{0}.
#' @param icp.batch.size A positive integer that specifies how many cells 
#' to randomly select. It behaves differently depending on \code{build.train.set}. 
#' If \code{build.train.set=FALSE}, it randomly samples cells for each ICP run 
#' from the complete dataset. If \code{build.train.set=TRUE}, it randomly samples 
#' cells once, before building the training set with the sampled cells (per batch 
#' if \code{batch.label} different than \code{NULL}). Default is \code{Inf}, 
#' which means using all cells.
#' @param train.with.bnn Train data with batch nearest neighbors. Default is 
#' \code{TRUE}. Only used if \code{batch.label} is given.   
#' @param train.k.nn Train data with batch nearest neighbors using \code{k} 
#' nearest neighbors. Default is \code{10}. Only used if \code{train.with.bnn} 
#' is \code{TRUE} and \code{train.k.nn.prop} is \code{NULL}.
#' @param train.k.nn.prop A numeric (higher than 0 and lower than 1) corresponding 
#' to the fraction of cells per cluster to use as \code{train.k.nn} nearest 
#' neighbors. If \code{NULL} the number of \code{train.k.nn} nearest neighbors 
#' is equal to \code{train.k.nn}. If given, \code{train.k.nn} parameter is ignored 
#' and \code{train.k.nn} is calculated based on \code{train.k.nn.prop}. By default 
#' \code{0.3} meaning that 30% of the cells are used. A vector with different 
#' proportions for the different divisive clustering rounds can be given, otherwise 
#' the same value is given for all.    
#' @param build.train.set Logical specifying if a training set should be built 
#' from the data or the whole data should be used for training. By default 
#' \code{TRUE}.
#' @param build.train.params A list of parameters to be passed to the function
#' \code{AggregateDataByBatch()}. Only provided if \code{build.train.set} is \code{TRUE}. 
#' @param scale.by A character specifying if the data should be scaled by \code{cell}
#' or by \code{feature} before training. Default is \code{NULL}, i.e., the data is 
#' not scaled before training.
#' @param use.cluster.seed Should the same starting clustering result be provided
#' to ensure more reproducible results (logical). If \code{FALSE}, each ICP run 
#' starts with a total random clustering and, thus, independent clustering. By 
#' default \code{TRUE}, i.e., the same clustering result is provided based on PCA 
#' density sampling. If \code{batch.label} different than \code{NULL}, the PCA
#' density sampling is performed in a batch wise manner.  
#' @param divisive.method Divisive method (character). One of \code{"random"} 
#' (randomly sample two clusters out of every cluster previously found),
#' \code{"cluster"} or \code{"cluster.batch"} (sample two clusters out of every 
#' cluster previously found based on the cluster probability distribution across
#' batches or per batch). By default \code{"cluster.batch"}. If \code{batch.label} 
#' is \code{NULL}, it is automatically set to \code{cluster}. It can be set to 
#' \code{random} if explicitly provided. 
#' @param allow.free.k Allow free \code{k} (logical). Allow ICP algorithm to 
#' decrease the \code{k} given in case it does not find \code{k} target clusters. 
#' By default \code{TRUE}. 
#' @param ari.cutoff Include ICP models and probability tables with an Adjusted 
#' Rand Index higher than \code{ari.cutoff} (numeric). By default \code{0.3}. A
#' value that can range between 0 (include all) and lower than 1.   
#' @param verbose A logical value to print verbose during the ICP run in case 
#' of parallelization, i.e., 'threads' different than \code{1}. Default 'FALSE'. 
#'
#' @name RunParallelDivisiveICP
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @keywords iterative clustering projection ICP logistic regression LIBLINEAR
#'
#' @importFrom parallel makeCluster detectCores stopCluster clusterCall
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doSNOW registerDoSNOW
#' @import Matrix
#' @import aricode
#' @import LiblineaR
#' @import SparseM
#' @importFrom SingleCellExperiment logcounts
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
#'                               threads = 2)
#' 
#' # Integrated PCA
#' set.seed(125) # to ensure reproducibility for the default 'irlba' method
#' sce <- RunPCA(object = sce, assay.name = "joint.probability", p = 10)
#' 
#' # Plot result 
#' cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
#'                               legend.nrow = 1),
#'                    PlotDimRed(object = sce, color.by = "Species", 
#'                               legend.nrow = 1), ncol = 2)
#'                    
RunParallelDivisiveICP.SingleCellExperiment <- function(object, batch.label, 
                                                        k, d, L, r, C,
                                                        reg.type, max.iter,
                                                        threads, icp.batch.size, 
                                                        train.with.bnn, train.k.nn,
                                                        train.k.nn.prop,
                                                        build.train.set, 
                                                        build.train.params,
                                                        scale.by, use.cluster.seed,
                                                        divisive.method,
                                                        allow.free.k,
                                                        ari.cutoff,
                                                        verbose) {
    
    if (!is(object,"SingleCellExperiment")) {
        stop("object must be of 'sce' class")
        return(object)
    }
    
    if (!is.null(batch.label)) {
        metadata(object)$coralysis$batch.label <- batch.label
    } else {
        if (divisive.method == "cluster.batch") {
            cat("WARNING: Setting 'divisive.method' to 'cluster' as 'batch.label=NULL'.", 
                "\nIf 'batch.label=NULL', 'divisive.method' can be one of: 'cluster', 'random'.", 
                "\n")
            divisive.method <- "cluster"
        }
    }
    
    if (!(is.numeric(k) && (k > 1) && ((log2(k) %% 1) == 0))) {
        stop("k must be a positive integer power of two")
    } else {
        metadata(object)$coralysis$k <- k
    }
    
    if (!is.numeric(d) | d >= 1 | d <= 0) {
        stop("d must be a numeric and in the range of (0,1)")
    } else {
        metadata(object)$coralysis$d <- d
    }
    
    if (!is.numeric(L) | L <= 0 | L%%1!=0) {
        stop("L must be a positive integer and greater than 0")
    } else {
        metadata(object)$coralysis$L <- L
    }
    
    if (!is.numeric(r) | r <= 0 | r%%1!=0) {
        stop("r must be a positive integer and greater than 0")
    } else {
        metadata(object)$coralysis$r <- r
    }
    
    if (!is.numeric(C) | C <= 0) {
        stop("C must be a numeric and greater than 0")
    } else {
        metadata(object)$coralysis$C <- C
    }
    
    if (!is.character(reg.type) | (reg.type != "L1" & reg.type != "L2")) {
        stop("reg.type parameter must be either 'L1' or 'L2'")
    } else {
        metadata(object)$coralysis$reg.type <- reg.type
    }
    
    if (!is.numeric(max.iter) | max.iter <= 0 | max.iter%%1 != 0) {
        stop("max.iter must be a positive integer and greater than 0")
    } else {
        metadata(object)$coralysis$max.iter <- max.iter
    }
    
    if (!is.numeric(threads) | threads < 0 | threads%%1 != 0) {
        stop("threads must be a positive integer or 0 (0 = use all available - 1)")
    } else {
        metadata(object)$coralysis$threads <- threads
    }
    
    if (!is.infinite(icp.batch.size)) {
        if (!is.numeric(icp.batch.size) | icp.batch.size <= 2 | icp.batch.size%%1 != 0) {
            stop("icp.batch.size must be a positive integer > 2 or Inf (0 = use all cells in ICP)")
        }
    } else {
        metadata(object)$coralysis$icp.batch.size <- icp.batch.size
    }
    
    if (!is.logical(build.train.set)) {
        stop("'build.train.set' must be logical")
    } else {
        metadata(object)$coralysis$build.train.set <- build.train.set
        
        if (build.train.set) {
            if (!is.list(build.train.params)) {
                stop("'build.train.params' must be a list")
            }
            metadata(object)$coralysis$build.train.params <- build.train.params
        }
    }
    
    if(!(is.null(scale.by) || ((length(scale.by)==1) && (scale.by %in% c("cell", "feature"))))) {
        stop("'scale.by' must be NULL or one of 'cell'/'feature'")
    } else {
        metadata(object)$coralysis$scale.by <- scale.by
    }
    
    if (!is.logical(use.cluster.seed)) {
        stop("'use.cluster.seed' must be logical")
    } else {
        metadata(object)$coralysis$use.cluster.seed <- use.cluster.seed
    }
    
    if (!is.logical(allow.free.k)) {
        stop("'allow.free.k' must be logical")
    } else {
        metadata(object)$coralysis$allow.free.k <- allow.free.k
    }
    
    if (!(is.numeric(ari.cutoff) && (length(ari.cutoff)==1) && ((ari.cutoff >= 0) || (ari.cutoff < 1)))) {
        stop("'ari.cutoff' must be a numeric value ranging between [0,1[")
    } else {
        metadata(object)$coralysis$ari.cutoff <- ari.cutoff
    }
    
    if (!(is.character(divisive.method) && (length(divisive.method)==1) && (divisive.method %in% c("random", "cluster", "cluster.batch")))) {
        stop("'divisive.method' must be a character among one of 'random'/'cluster'/'cluster.batch'")
    } else {
        metadata(object)$coralysis$divisive.method <- divisive.method
    }
    
    if (!is.logical(train.with.bnn)) {
        stop("'train.with.bnn' must be logical")
    } else {
        metadata(object)$coralysis$train.with.bnn <- train.with.bnn
        if (!(is.null(train.k.nn.prop) || (is.numeric(train.k.nn.prop) && (length(train.k.nn.prop)==1) && ((train.k.nn.prop > 0) || (train.k.nn.prop < 1))))) {
            stop("'train.k.nn.prop' must be NULL or ranging between ]0,1[")
        } else {
            metadata(object)$coralysis$train.k.nn.prop <- train.k.nn.prop
            if (is.null(train.k.nn.prop) && !(is.numeric(train.k.nn) && (length(train.k.nn)==1) && (train.k.nn%%1 == 0))) {
                stop("'train.k.nn' must be an integer")
            } else {
                metadata(object)$coralysis$train.k.nn <- train.k.nn
            }
        }
    }
    metadata(object)$coralysis$divisive.icp <- TRUE
    
    if (build.train.set) {
        if (!is.infinite(icp.batch.size)) { # check if dataset should be sampled
            if (!is.null(batch.label)) {
                # sample cells per batch dataset
                message("\nDownsampling dataset by batch using random sampling.")
                cellidx.by.batch <- split(x = 1:ncol(object), f = object[[batch.label]])
                ncells.by.batch <- unlist(lapply(X = cellidx.by.batch, FUN = length))
                ncells.by.batch[ncells.by.batch > icp.batch.size] <- icp.batch.size
                if (any(ncells.by.batch < icp.batch.size)) { # check if any batch dataset has lower number of cells than icp.batch.size
                    cat("WARNING: the number of cells in one or more batch dataset(s)", 
                        "\nis lower than the 'icp.batch.size':", icp.batch.size,
                        "\nUsing all the available cells instead:\n", 
                        ncells.by.batch[ncells.by.batch < icp.batch.size])
                } 
                pick.cells <- lapply(X = seq_along(ncells.by.batch), FUN = function(x) { # sample cell barcode names per batch dataset
                    sample(x = cellidx.by.batch[[x]], size = ncells.by.batch[[x]], replace = FALSE) 
                })
                pick.cells <- unlist(pick.cells)
            } else {
                # Sample complete dataset
                message("\nDownsampling overall dataset using random sampling.")
                if (ncol(object) < icp.batch.size) {
                    message(cat("WARNING: the number of cells in current batch is", 
                                ncol(object), "lower than the 'icp.batch.size' -", 
                                icp.batch.size, "\nUsing all the available cells instead:", 
                                ncol(object)))
                    icp.batch.size <- ncol(object)
                }
                pick.cells <- sample(x = 1:ncol(object), size = icp.batch.size, replace = FALSE) 
            }
            icp.batch.size <- Inf # prevent recurrent sampling below in 'RunDivisiveICP'
        } else { # use all the cells from all the batch datasets
            pick.cells <- 1:ncol(object)
        }
        message("\nBuilding training set...")
        build.train.params <- c(list(object = object[,pick.cells], batch.label = batch.label), build.train.params)
        clustered.object <- do.call(AggregateDataByBatch, build.train.params)
        message("Training set successfully built.")
        dataset <- t(logcounts(clustered.object))
    } else {
        dataset <- t(logcounts(object))
    }
    
    batch.name <- batch.label
    if (!is.null(batch.label)) {
        if (build.train.set) {
            batch.label <- as.character(clustered.object[["batch"]])
            names(batch.label) <- colnames(clustered.object)
        } else {
            batch.label <- as.character(object[[batch.label]])
            names(batch.label) <- colnames(object)
        }
    }
    
    if (use.cluster.seed) {
        message("\nComputing cluster seed.")
        cluster.seed <- SamplePCACells(data = dataset, batch = batch.label, 
                                       q.split = 0.5, p=30, use.pc="PC1", 
                                       center=TRUE, scale.=TRUE)
    } else {
        cluster.seed <- NULL
    }
    metadata(object)$coralysis$scale.by <- scale.by
    if (!is.null(scale.by)) {
        message(cat("Scaling data by", scale.by, ".\n"))
        if (scale.by=="cell") {
            dataset <- Scale(x = as(dataset, "sparseMatrix"), scale.by="row")
        }
        if (scale.by=="feature") {
            dataset <- ScaleByBatch(x = dataset, batch = batch.label)
        }
    }
    parallelism <- TRUE

    if (threads == 0) {
        cl <- makeCluster(detectCores(logical=TRUE)-1)
        # registerDoParallel(cl)
        registerDoSNOW(cl)
    } else if(threads == 1) {
        message("Parallelism disabled, because threads = 1")
        parallelism <- FALSE
    } else {
        if (verbose) {
            cl <- makeCluster(threads, outfile="")
        } else {
            cl <- makeCluster(threads) 
        }
        # registerDoParallel(cl)
        registerDoSNOW(cl)
        clusterCall(cl, function(x) .libPaths(x), .libPaths()) # Exporting .libPaths from master to the workers
    }
    
    message("\nInitializing divisive ICP clustering...\n")
    if (parallelism) {
        pb <- txtProgressBar(min = 1, max = L, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        out <- foreach(task = seq_len(L),
                       .verbose = FALSE,
                       .combine = list,
                       .maxcombine = 1000,
                       .inorder = TRUE,
                       .multicombine = TRUE,
                       .packages=c("Coralysis", "parallel"),
                       .noexport = "task",
                       .options.snow = opts)  %dorng% {
                           tryCatch(expr = {
                               task <- task
                               message(paste0("\nICP run: ", task))
                               RunDivisiveICP(normalized.data = dataset, batch.label = batch.label, 
                                              k = k, d = d, r = r, C = C, reg.type = reg.type, 
                                              max.iter = max.iter, icp.batch.size = icp.batch.size, 
                                              train.with.bnn = train.with.bnn, train.k.nn = train.k.nn, 
                                              train.k.nn.prop = train.k.nn.prop, cluster.seed = cluster.seed, 
                                              divisive.method = divisive.method, allow.free.k = allow.free.k, 
                                              ari.cutoff = ari.cutoff)
                           }, error = function(e){ # Stop progress bar & workers if 'foreach()' loop terminates/exit with error
                               message("'foreach()' loop terminated unexpectedly.\nPlease read the error message or use the 'verbose=TRUE' option.\nShutting down workers...")
                               close(pb)
                               stopCluster(cl)
                           })
                       }
        close(pb)
        # stop local cluster
        stopCluster(cl)
    } else {
        out <- list()
        for (l in seq_len(L)) {
            try({
                message(paste0("ICP run: ",l))
                res <- RunDivisiveICP(normalized.data = dataset, batch.label = batch.label, 
                                      k = k, d = d, r = r, C = C, reg.type = reg.type, 
                                      max.iter = max.iter, icp.batch.size = icp.batch.size, 
                                      train.with.bnn = train.with.bnn, train.k.nn = train.k.nn, 
                                      train.k.nn.prop = train.k.nn.prop, cluster.seed = cluster.seed, 
                                      divisive.method = divisive.method, allow.free.k = allow.free.k, 
                                      ari.cutoff = ari.cutoff)
                out[[l]] <- res
            })
        }
    }
    message("\nDivisive ICP clustering completed successfully.")
    metadata(object)$coralysis$metrics <-
        lapply(out, function(x) x$metrics)
    metadata(object)$coralysis$metrics <- unlist(metadata(object)$coralysis$metrics, recursive = FALSE)
    metadata(object)$coralysis$models <-
        lapply(out, function(x) x$model)
    metadata(object)$coralysis$models <- unlist(metadata(object)$coralysis$models, recursive = FALSE)
    
    message("\nPredicting cell cluster probabilities using ICP models...")
    if (build.train.set) {
        test.data <- t(logcounts(object))
        colnames(test.data) <- paste0("W", 1:ncol(test.data))
        if (!is.null(scale.by)) {
            if (scale.by=="cell") {
                test.data <- Scale(x = as(test.data, "sparseMatrix"), scale.by="row")
            }
            if (scale.by=="feature") {
                batch.label <- as.character(object[[batch.name]]) 
                names(batch.label) <- colnames(object)
                test.data <- ScaleByBatch(x = test.data, batch = batch.label)
            }
        }
        metadata(object)$coralysis$joint.probability <- 
            lapply(metadata(object)$coralysis$models, function(x) {
                predict(x, test.data, proba=TRUE)$probabilities
            })
    } else {
        metadata(object)$coralysis$joint.probability <-
            unlist(lapply(out, function(x) x$probabilities), recursive = FALSE)
    }
    message("Prediction of cell cluster probabilities completed successfully.")
    message("\nMulti-level integration completed successfully.")
    
    return(object)
}
#' @rdname RunParallelDivisiveICP
#' @aliases RunParallelDivisiveICP
setMethod("RunParallelDivisiveICP", signature(object = "SingleCellExperiment"),
          RunParallelDivisiveICP.SingleCellExperiment)


#' @title Aggregates feature expression by cell clusters, per batch if provided.
#'
#' @description The function aggregates feature expression by cell clusters, per 
#' batch if provided.
#'
#' @param object An object of \code{SingleCellExperiment} class.
#' @param batch.label A variable name (of class \code{character}) available 
#' in the cell metadata \code{colData(object)} with the batch labels (\code{character} 
#' or \code{factor}) to use. The variable provided must not contain \code{NAs}.
#' By default \code{NULL}, i.e., cells are sampled evenly regardless their batch. 
#' @param batch.label Cluster identities vector corresponding to the cells in 
#' \code{mtx}.
#' @param nhvg Integer of the number of highly variable features to select. By default 
#' \code{2000}. 
#' @param p Integer. By default \code{30}. 
#' @param ... Parameters to be passed to \code{ClusterCells()} function. 
#'
#' @name AggregateDataByBatch
#' 
#' @return A SingleCellExperiment object with feature expression aggregated by clusters.
#'
#' @keywords aggregated feature expression batches
#'
#' @importFrom SingleCellExperiment logcounts cbind
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom scran getTopHVGs
#' @importFrom irlba prcomp_irlba
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
#' # Run with a batch
#' set.seed(1204)
#' sce <- AggregateDataByBatch(object = sce, batch.label = "batch")
#' logcounts(sce)[1:10,1:10]
#' head(metadata(sce)$clusters)
#' 
#' # Run without a batch
#' set.seed(1204)
#' sce <- AggregateDataByBatch(object = sce, batch.label = NULL)
#' logcounts(sce)[1:10,1:10]
#' head(metadata(sce)$clusters)
#' } 
#' 
AggregateDataByBatch.SingleCellExperiment <- function(object, batch.label, 
                                                      nhvg, p, ...) {
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), 
              (is.null(batch.label) || (is.character(batch.label) && (batch.label %in% colnames(colData(object))))), 
              (is.numeric(nhvg) && (nhvg%%1 == 0) && (nhvg <= nrow(object))), 
              (is.numeric(p) && (p%%1 == 0) && (p <= nrow(object))))
    
    if (is.null(batch.label)) {
        batch.label <- "batch"
        object[[batch.label]] <- "single"
    }
    batch <- as.character(colData(object)[[batch.label]])
    batch.names <- unique(batch)
    names(batch.names) <- batch.names
    sce.batch <- lapply(X = batch.names, FUN = function(b) {
        object[,batch==b]
    })
    top.hvg <- lapply(X = sce.batch, FUN = function(b) {
        getTopHVGs(b, n=nhvg)
    })
    pca.batch <- lapply(X = batch.names, FUN = function(b) {
        prcomp_irlba(t(logcounts(sce.batch[[b]][top.hvg[[b]],])),
                     scale.=TRUE, center=TRUE, n=p)
    })
    meta.data <- data.frame()
    sce.batch.clusters <- list()
    for (b in batch.names) {
        reducedDim(x=sce.batch[[b]], type="PCA") <- pca.batch[[b]]$x
        sce.batch[[b]] <- ClusterCells(object = sce.batch[[b]], ...)
        clusters.mean <- AggregateClusterExpression(mtx = logcounts(sce.batch[[b]]),
                                                    cluster=as.character(metadata(sce.batch[[b]])$clusters@cluster))
        colnames(clusters.mean) <- paste(colnames(clusters.mean), b, sep="_")
        sce.batch.clusters[[b]] <- SingleCellExperiment(assays=list(logcounts=clusters.mean), 
                                                        colData=DataFrame(batch=rep(b, ncol(clusters.mean)), 
                                                                          row.names = colnames(clusters.mean)))
        tmp.meta.data <- data.frame("cluster" = metadata(sce.batch[[b]])$clusters@cluster, "batch" = b)
        meta.data <- rbind(meta.data, tmp.meta.data)
    }
    sce <- do.call(cbind, sce.batch.clusters)
    metadata(sce)$clusters <- meta.data[colnames(object),]
    return(sce)
}
#' @rdname AggregateDataByBatch
#' @aliases AggregateDataByBatch
setMethod("AggregateDataByBatch", signature(object = "SingleCellExperiment"),
          AggregateDataByBatch.SingleCellExperiment)
