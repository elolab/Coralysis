#' @title Run ICP divisive clustering parallerly with increasing number of Ks
#'
#' @description
#' This functions runs in parallel \code{L} ICP runs, which is the computational
#' bottleneck of ILoReg2. With ~ 3,000 cells this step should be completed
#' in ~ 2 h and ~ 1 h with 3 and 12 logical processors (threads), respectively.
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
#' Decreasing leads to more stringent feature selection, i.e. less genes are
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
#' if \code{batch.label=NULL}).
#' Default is \code{Inf}, which means using all cells.
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
#' \code{AggregateDataByBatch()}.
#' @param scale.by A character specifying if the data should be scaled by \code{cell}
#' or by \code{gene} before training. Default is \code{NULL}, i.e., the data is 
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
#' @return an object of \code{SingleCellExperiment} class
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
        stop("object must of 'sce' class")
        return(object)
    }
    
    if (!is.null(batch.label)) {
        metadata(object)$iloreg$batch.label <- batch.label
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
        metadata(object)$iloreg$k <- k
    }
    
    if (!is.numeric(d) | d >= 1 | d <= 0)
    {
        stop("d must be a numeric and in the range of (0,1)")
    } else {
        metadata(object)$iloreg$d <- d
    }
    
    if (!is.numeric(L) | L <= 0 | L%%1!=0)
    {
        stop("L must be a positive integer and greater than 0")
    } else {
        metadata(object)$iloreg$L <- L
    }
    
    if (!is.numeric(r) | r <= 0 | r%%1!=0)
    {
        stop("r must be a positive integer and greater than 0")
    } else {
        metadata(object)$iloreg$r <- r
    }
    
    if (!is.numeric(C) | C <= 0)
    {
        stop("C must be a numeric and greater than 0")
    } else {
        metadata(object)$iloreg$C <- C
    }
    
    if (!is.character(reg.type) | (reg.type != "L1" & reg.type != "L2"))
    {
        stop("reg.type parameter must be either 'L1' or 'L2'")
    } else {
        metadata(object)$iloreg$reg.type <- reg.type
    }
    
    if (!is.numeric(max.iter) | max.iter <= 0 | max.iter%%1 != 0)
    {
        stop("max.iter must be a positive integer and greater than 0")
    } else {
        metadata(object)$iloreg$max.iter <- max.iter
    }
    
    if (!is.numeric(threads) | threads < 0 | threads%%1 != 0)
    {
        stop("threads must be a positive integer or 0 (0 = use all available - 1)")
    } else {
        metadata(object)$iloreg$threads <- threads
    }
    
    if (!is.infinite(icp.batch.size))
    {
        if (!is.numeric(icp.batch.size) | icp.batch.size <= 2 | icp.batch.size%%1 != 0)
        {
            stop("icp.batch.size must be a positive integer > 2 or Inf (0 = use all cells in ICP)")
        } else {
            metadata(object)$iloreg$icp.batch.size <- icp.batch.size
        }
    }
    metadata(object)$iloreg$divisive.icp <- TRUE
    
    if (build.train.set) {
        if (!is.infinite(icp.batch.size)) { # check if dataset should be sampled
            if (!is.null(batch.label)) {
                # sample cells per batch dataset
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
        build.train.params <- c(list(object = object[,pick.cells], batch.label = batch.label), build.train.params)
        clustered.object <- do.call(AggregateDataByBatch, build.train.params)
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
        cluster.seed <- SamplePCACells(data = dataset, batch = batch.label, 
                                       q.split = 0.5, p=30, use.pc="PC1", 
                                       center=TRUE, scale.=TRUE)
    } else {
        cluster.seed <- NULL
    }
    metadata(object)$iloreg$scale.by <- scale.by
    if (!is.null(scale.by)) {
        if (scale.by=="cell") {
            dataset <- Scale(x = as(dataset, "sparseMatrix"), scale.by="row")
        }
        if (scale.by=="gene"){
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
    
    if (parallelism) {
        pb <- txtProgressBar(min = 1, max = L, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        out <- foreach(task = seq_len(L),
                       .verbose = FALSE,
                       .combine = list,
                       .maxcombine = 1000,
                       .inorder = FALSE,
                       .multicombine = TRUE,
                       .packages=c("ILoReg2", "parallel"),
                       .options.snow = opts)  %dorng% {
                           tryCatch(expr = {
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
    metadata(object)$iloreg$metrics <-
        lapply(out, function(x) x$metrics)
    metadata(object)$iloreg$metrics <- unlist(metadata(object)$iloreg$metrics, recursive = FALSE)
    metadata(object)$iloreg$models <-
        lapply(out, function(x) x$model)
    metadata(object)$iloreg$models <- unlist(metadata(object)$iloreg$models, recursive = FALSE)
    
    if (build.train.set) {
        test.data <- t(logcounts(object))
        colnames(test.data) <- paste0("W", 1:ncol(test.data))
        if (!is.null(scale.by)) {
            if (scale.by=="cell") {
                test.data <- Scale(x = as(test.data, "sparseMatrix"), scale.by="row")
            }
            if (scale.by=="gene") {
                batch.label <- as.character(object[[batch.name]]) 
                names(batch.label) <- colnames(object)
                test.data <- ScaleByBatch(x = test.data, batch = batch.label)
            }
        }
        metadata(object)$iloreg$joint.probability <- 
            lapply(metadata(object)$iloreg$models, function(x) {
                predict(x, test.data, proba=TRUE)$probabilities
            })
    } else {
        metadata(object)$iloreg$joint.probability <-
            unlist(lapply(out, function(x) x$probabilities), recursive = FALSE)
    }

    return(object)
}
#' @rdname RunParallelDivisiveICP
#' @aliases RunParallelDivisiveICP
setMethod("RunParallelDivisiveICP", signature(object = "SingleCellExperiment"),
          RunParallelDivisiveICP.SingleCellExperiment)

#' @title Iterative Clustering Projection (ICP) divisive clustering
#'
#' @description
#' The function implements Iterative Clustering Projection (ICP): a
#' supervised learning-based clustering, which maximizes clustering similarity
#' between the clustering and its projection by logistic regression.
#'
#' @param normalized.data A sparse matrix (dgCMatrix) containing
#' normalized gene expression data with cells in rows and genes in columns.
#' Default is \code{NULL}.
#' @param batch.label A character vector with batch labels corresponding to the cells
#' given in \code{normalized.data}. The character batch labels need to be named
#' with the cells names given in the rows of \code{normalized.data}. 
#' By default \code{NULL}, i.e., cells are sampled evenly regardless their batch. 
#' @param k A positive integer power of two, i.e., \code{2**n}, where \code{n>0},
#' specifying the number of clusters in the last Iterative Clustering Projection (ICP)
#' round. Decreasing \code{k} leads to smaller cell populations diversity and vice versa. 
#' Default is \code{16}, i.e., the divisive clustering 2 -> 4 -> 8 -> 16 is performed. 
#' @param d A numeric that defines how many cells per cluster should be
#' down- and oversampled (d in ceiling(N/k*d)), when stratified.downsampling=FALSE,
#' or what fraction should be downsampled in the stratified approach
#' ,stratified.downsampling=TRUE. Default is \code{0.3}.
#' @param r A positive integer that denotes the number of reiterations
#' performed until the algorithm stops. Default is \code{5}.
#' @param C Cost of constraints violation (\code{C}) for L1-regulatization.
#' Default is \code{0.3}.
#' @param reg.type "L1" for LASSO and "L2" for Ridge. Default is "L1".
#' @param max.iter A positive integer that denotes the maximum number of
#' iterations performed until the algorithm ends. Default is \code{200}.
#' @param icp.batch.size A positive integer that specifies how many cells 
#' to randomly select for each ICP run from the complete data set. 
#' This is a new feature intended to speed up the process
#' with larger data sets. Default is \code{Inf}, which means using all cells.
#' @param train.with.bnn Train data with batch nearest neighbors. Default is 
#' \code{TRUE}. Only used if \code{batch.label} is given.   
#' @param train.k.nn Train data with batch nearest neighbors using \code{k} 
#' nearest neighbors. Default is \code{10}. Only used if \code{train.with.bnn} 
#' is \code{TRUE}.  
#' @param train.k.nn.prop A numeric (higher than 0 and lower than 1) corresponding 
#' to the fraction of cells per cluster to use as \code{train.k.nn} nearest 
#' neighbors. Default is \code{NULL} meaning that the number of \code{train.k.nn} 
#' nearest neighbors is equal to \code{train.k.nn}. If given, \code{train.k.nn} 
#' parameter is ignored and \code{train.k.nn} is calculated based on 
#' \code{train.k.nn.prop}. A vector with different proportions for the different
#' divisive clustering rounds can be given, otherwise the same value is given for 
#' all.   
#' @param cluster.seed A cluster seed to start and guide the clustering to more 
#' reproducible clusterings across runs (factor). Default is \code{NULL}. Otherwise, 
#' a random clustering takes place to start divisive clustering with ICP. 
#' @param divisive.method Divisive method (character). One of \code{"random"} 
#' (randomly sample two clusters out of every cluster previously found),
#' \code{"cluster"} or \code{"cluster.batch"} (sample two clusters out of every 
#' cluster previously found based on the cluster probability distribution across
#' batches or per batch). By default \code{"random"}. 
#' @param allow.free.k Allow free \code{k} (logical). Allow ICP algorithm to 
#' decrease the \code{k} given in case it does not find \code{k} target clusters. 
#' By default \code{FALSE}. 
#' @param ari.cutoff Include ICP models and probability tables with an Adjusted 
#' Rand Index higher than \code{ari.cutoff} (numeric). By default \code{0.5}. A
#' value that can range between 0 (include all) and lower than 1.   
#' 
#' @return A list that includes the probability matrix and the clustering
#' similarity measures: ARI, NMI, etc.
#'
#' @keywords iterative clustering projection ICP clustering
#'
#' @import Matrix
#' @importFrom aricode clustComp
#' @import LiblineaR
#' @import SparseM
#'
RunDivisiveICP <- function(normalized.data = NULL, batch.label = NULL, 
                           k = 16, d = 0.3, r = 5, C = 5,
                           reg.type = "L1", max.iter = 200, 
                           icp.batch.size=Inf, train.with.bnn = TRUE, 
                           train.k.nn = 10, train.k.nn.prop = NULL, 
                           cluster.seed = NULL, divisive.method = "random", 
                           allow.free.k = FALSE, ari.cutoff = 0.5) {
    
    metrics <- NULL
    idents <- list()
    iterations <- 1
    
    if ((!is.infinite(icp.batch.size)) & icp.batch.size > 2) {
        if (nrow(normalized.data) < icp.batch.size) {
            message(cat("WARNING: the number of cells in current batch is", 
                        nrow(normalized.data), "lower than the 'icp.batch.size' -", 
                        icp.batch.size, "\nUsing all the available cells instead:", 
                        nrow(normalized.data)))
            icp.batch.size <- nrow(normalized.data)
        }
        normalized_data_whole <- normalized.data
        randinds_batch <- sample(seq_len(nrow(normalized.data)),
                                 size = icp.batch.size,replace = FALSE)
        normalized.data <- normalized.data[randinds_batch,]
    }
    
    probs <- res_model <- res_metrics <- preds <- list()
    i <- 0
    Ks <- 2^seq(from=1, to=log2(k), by=1)
    if (length(train.k.nn.prop)==1) {
        train.k.nn.prop <- rep(train.k.nn.prop, length(Ks))
    }
    for (k in Ks) {
        first_round <- TRUE
        i <- i + 1
        while (TRUE) {
            
            # Step 1: initialize clustering (ident_1) randomly, ARI=0 and r=0
            if (first_round) {
                if (k == 2) {
                    if (!is.null(cluster.seed) & is.factor(cluster.seed)) {
                        ident_1 <- cluster.seed
                    } else {
                        ident_1 <- factor(sample(seq_len(k), nrow(normalized.data), replace = TRUE))
                    }
                } else {
                    if (divisive.method=="random") {
                        ident_1 <- RandomlyDivisiveClustering(cluster=preds[[i-1]], k=2)
                    } 
                    if (divisive.method=="cluster") {
                        ident_1 <- SampleClusterProbs(cluster = preds[[i-1]], probs = probs[[i-1]], q.split = 0.5)
                    }
                    if (divisive.method=="cluster.batch") {
                        ident_1 <- SampleClusterBatchProbs(cluster = preds[[i-1]], probs = probs[[i-1]], 
                                                           batch = batch.label, q.split = 0.5)
                    }
                }
                names(ident_1) <- row.names(normalized.data)
                idents[[1]] <- ident_1
                ari <- 0
                reiterations <- 0
            }
            
            # Step 2: train logistic regression model
            if (!is.null(batch.label) & train.with.bnn & !first_round) {
                training_ident_subset <- unlist(
                    FindClusterBatchKNN(preds = res_prediction$predictions, 
                                        probs = res_prediction$probabilities,
                                        batch = batch.label, 
                                        k = train.k.nn, 
                                        k.prop = train.k.nn.prop[i])
                )
            } else {
                training_ident_subset <- NULL
            }
            res <- LogisticRegression(training.sparse.matrix = normalized.data,
                                      training.ident = ident_1, C = C,
                                      reg.type=reg.type,
                                      test.sparse.matrix = normalized.data, d=d, 
                                      batch.label = batch.label, 
                                      training_ident_subset = training_ident_subset)
            
            res_prediction <- res$prediction
            
            names(res_prediction$predictions) <- row.names(normalized.data)
            rownames(res_prediction$probabilities) <- row.names(normalized.data)
            
            message(paste0("probability matrix dimensions = ",paste(dim(res_prediction$probabilities),collapse = " ")))
            
            # Projected clusters
            ident_2 <- res_prediction$predictions
            
            # Safety procedure: If k drops to 1, start from the beginning.
            # However, k should NOT decrease during the iteration when
            # the down- and oversampling approach is used for the balancing training data.
            if (nlevels(factor(as.character(ident_2))) < k) {
                if (allow.free.k & (nlevels(factor(as.character(ident_2)))>3)) {
                    message(paste("k", k, "decreased to", nlevels(factor(as.character(ident_2)))))
                    k <- nlevels(factor(as.character(ident_2)))
                } else {
                    message(paste0("k decreased, starting from the beginning... ",
                                   "consider increasing d to 0.5 and C to 1 ",
                                   "or increasing the ICP batch size ",
                                   "and check the input data ",
                                   "(scaled dense data might cause problems)"))
                    first_round <- TRUE
                    metrics <- NULL
                    idents <- list()
                    iterations <- 1
                    next
                }
            }
            
            # Step 3: compare clustering similarity between clustering and projection
            message(paste0("EPOCH: ",iterations))
            message(paste0("current clustering = ",paste(table(ident_1),collapse = " ")))
            message(paste0("projected clustering = ",paste(table(ident_2),collapse = " ")))
            
            comp_clust <- clustComp(c1 = ident_1, c2 = ident_2)
            
            if(first_round & comp_clust$ARI <= 0) {
                next
            }
            
            # Step 3.1: If ARI did not increase, reiterate
            if (comp_clust$ARI <= ari & !(first_round)) {
                reiterations <- reiterations + 1
            } else { # Step 3.2: If ARI increased, proceed to next iteration round
                # Update clustering to the predicted clusters
                message(paste0("ARI=",as.character(comp_clust$ARI)))
                ident_1 <- ident_2
                first_round <- FALSE
                metrics <- cbind(metrics, comp_clust)
                iterations = iterations + 1
                idents[[iterations]] <- ident_2
                ari <- comp_clust$ARI
                reiterations <- 0
                # save model & cluster probs
                res_model[[i]] <- res$model
                probs[[i]] <- res_prediction$probabilities
                preds[[i]] <- ident_2
            }
            # Step 4: If the maximum number of reiterations or iterations was
            # reached, break the while loop
            if (reiterations == r | iterations == max.iter){
                if (ari>ari.cutoff) {
                    break
                } else {
                    first_round <- TRUE
                    metrics <- NULL
                    idents <- list()
                    iterations <- 1
                    next
                }
            } 
        }
        res_metrics[[i]] <- metrics
        iterations <- 1
    }
    
    message("ICP converged at EPOCH ", iterations, ".\nMaximum ARI reached: ", as.character(ari),".")
    if (is.infinite(icp.batch.size)) {
        # Step 5: Return result
        return(list(probabilities=probs, metrics=res_metrics,model=res_model))
    } else {
        cat("projecting the whole data set...")
        res <- LogisticRegression(training.sparse.matrix = normalized.data,
                                  training.ident = ident_1, C = C,
                                  reg.type=reg.type,
                                  test.sparse.matrix = normalized_data_whole, d=d, 
                                  batch.label = batch.label, 
                                  training_ident_subset = training_ident_subset)
        
        res_prediction <- res$prediction
        res_model <- res$model
        
        names(res_prediction$predictions) <- row.names(normalized_data_whole)
        rownames(res_prediction$probabilities) <- row.names(normalized_data_whole)
        
        # Projected cluster probabilities
        probs <- res_prediction$probabilities
        
        message(" success!")
        
        message(paste(dim(probs),collapse = " "))
        
        return(list(probabilities=probs, metrics=metrics,model=res_model))
    }
}

#' @title Split randomly every cluster into k clusters
#'
#' @description
#' Splits randomly every cluster given into k clusters
#'
#' @param cluster A clustering result in which each cluster will be divided 
#' randomly into \code{k} clusters. 
#' @param k The number of clusters that each cluster given in \code{cluster} 
#' should be randomly divided into. 
#' @param cluster.names Names to name the \code{cluster} result given. By 
#' default \code{NULL}.
#' 
#' @return A clustering result where every cluster given was split randomly 
#' into \code{k} clusters. 
#'
#' @keywords random clustering
#'
RandomlyDivisiveClustering <- function(cluster, k, cluster.names = NULL) {
    clt.len <- 1:length(cluster) 
    clt.idx <- split(x = clt.len, f = cluster)
    clt.list <- lapply(X= clt.idx, function(x) {
        factor(sample(seq_len(k), length(x), replace = TRUE))
    })
    no.clts <- length(clt.list)
    clt.seq <- seq(from=k, to=k*no.clts, by=k)
    ii <- 0
    for (i in seq_along(clt.list)) {
        if (i > 1) {
            levels(clt.list[[i]]) <- (clt.seq[ii]+1):clt.seq[i]  
        }
        ii <- i
    }
    rand.divisive.clts <- unlist(clt.list)
    rand.divisive.clts <- rand.divisive.clts[order(unlist(clt.idx))]
    if (!is.null(cluster.names)) names(rand.divisive.clts) <- cluster.names
    return(rand.divisive.clts)
}

#' @title Sample cells based on principal components distribution
#'
#' @description
#' Samples cells based on their distributions along one principal component
#'
#' @param data Data to compute PCA and sample cells from. Rows and columns 
#' should represent cells and features/genes, respectively.   
#' @param batch Batch cell label identity (character) matching cells giving in 
#' \code{data}. Use \code{NULL} in the absence of batches. If the batch is given 
#' the cells are sampled in a batch wise manner, otherwise the cells are sampled 
#' without any grouping factor. By default is \code{NULL}.   
#' @param q.split Split (cell) batch principal component distribution by this 
#' quantile (numeric). By default {0.5}, i.e., median.  
#' @param p Number of principal components to compute (integer). By default 
#' \code{30}.
#' @param use.pc Which principal component should be used for sampling cells per
#' batch. By default \code{"PC1"}, i.e., first principal component is used. 
#' @param center Should the features/genes given in \code{data} centered before
#' performing the PCA (logical). By default \code{TRUE}. 
#' @param scale. Should the features/genes given in \code{data} scaled before
#' performing the PCA (logical). By default \code{TRUE}. 
#' 
#' @return A factor with cell cluster identities (two clusters). 
#'
#' @keywords sample PCA cells
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom stats quantile
#' @importFrom sparseMatrixStats colSds
#' 
SamplePCACells <- function(data, batch = NULL, q.split = 0.5, p=30, use.pc="PC1", 
                           center=TRUE, scale.=TRUE) {
    # Check batch
    if (is.null(batch)) {
        batch <- factor(rep("0", nrow(data)))
    }
    
    # Compute PCA
    genes.sd <- colSds(data) 
    genes.sd.diff0 <- which(genes.sd != 0)
    cell.names <- row.names(data)
    pca <- prcomp_irlba(x=data[,genes.sd.diff0], n=p, center = center, scale. = scale.)
    # Select PC of interest
    pc <- pca$x[,use.pc]
    names(pc) <- cell.names
    # Index by batch labels
    cell.batch <- split(x = 1:length(batch), f = batch)
    pc.batch <- split(x = pc, f = batch)
    # Split PC by quantile - median by default
    batch.names <- names(cell.batch)
    names(batch.names) <- batch.names
    q.batch <- lapply(X = pc.batch, FUN = function(x) quantile(x = x, probs = q.split))
    split.half.batch <- lapply(X = batch.names, FUN = function(x) {
        (pc.batch[[x]] <= q.batch[[x]])
    })
    first.half.batch <- lapply(X = batch.names, function(x) {
        names(pc.batch[[x]][split.half.batch[[x]]]) 
    })
    second.half.batch <- lapply(X = batch.names, function(x) {
        names(pc.batch[[x]][!split.half.batch[[x]]]) 
    })
    # Clustering 
    first.half.batch <- unlist(first.half.batch)
    second.half.batch <- unlist(second.half.batch)
    cluster <- factor(rep(1:2, c(length(first.half.batch), length(second.half.batch))))
    names(cluster) <- c(first.half.batch, second.half.batch)
    cluster <- factor(cluster[cell.names])
    return(cluster)
}

#' @title Sample cells based on cluster probabilities distribution
#'
#' @description
#' Samples cells based on cluster probabilities distribution
#'
#' @param cluster Clustering cell labels predicted by ICP (factor). 
#' @param probs Clustering probabilities predicted by ICP (matrix). 
#' @param q.split Split (cell) batch principal component distribution by this 
#' quantile (numeric). By default {0.5}, i.e., median.  

#' @return A factor with cell cluster identities. 
#'
#' @keywords sample cluster probabilities ICP
#' 
#' @importFrom stats quantile
#'
SampleClusterProbs <- function(cluster, probs, q.split = 0.5) {
    # Split cells & probs by cluster/prediction
    clts <- as.character(unique(cluster))
    names(clts) <- clts
    cell.idx <- 1:length(cluster)
    clt.cell.idx <- split(x = cell.idx, f = cluster)
    max.probs <- apply(X = probs, MARGIN = 1, FUN = function(x) max(x)) 
    clt.probs <- split(x = max.probs, f = cluster)
    # Quantile - median by default
    clt.q.split <- lapply(X = clt.probs, FUN = function(x) quantile(x, probs = q.split))
    clt.q.probs.split <- lapply(X = clts, FUN = function(x) (clt.probs[[x]] <= clt.q.split[[x]]))
    # Make new clusterings
    new.idx <- new.clt <- vector(mode="numeric", length = length(cluster))
    s.pos <- 1
    e.pos <- i <- 0
    for (clt in clts) {
        new.clt1 <- clt.cell.idx[[clt]][clt.q.probs.split[[clt]]]
        new.clt2 <- clt.cell.idx[[clt]][!clt.q.probs.split[[clt]]]
        make.clt <- rep(c(1+i, 2+i), c(length(new.clt1), length(new.clt2)))
        e.pos <- e.pos + length(make.clt) 
        new.idx[s.pos:e.pos] <- c(new.clt1, new.clt2)
        new.clt[s.pos:e.pos] <- make.clt
        i <- i + 2
        s.pos <- s.pos + length(make.clt) 
    }
    new.clt <- factor(new.clt[order(new.idx)])
    return(new.clt)
}

#' @title Sample cells based on cluster probabilities distribution batch wise
#'
#' @description
#' Samples cells based on cluster probabilities distribution batch wise
#'
#' @param cluster Clustering cell labels predicted by ICP (factor). 
#' @param probs Clustering probabilities predicted by ICP (matrix). 
#' @param batch Batch labels for the corresponding clusters (character or factor). 
#' @param q.split Split (cell) batch principal component distribution by this 
#' quantile (numeric). By default {0.5}, i.e., median.  

#' @return A factor with cell cluster identities. 
#'
#' @keywords sample cluster probabilities ICP
#' 
#' @importFrom stats quantile
#'
SampleClusterBatchProbs <- function(cluster, probs, batch, q.split = 0.5) {
    # Split cells & probs by cluster/prediction
    clts <- as.character(sort(unique(cluster)))
    names(clts) <- clts
    batches <- unique(as.character(batch))
    names(batches) <- batches
    cell.idx <- 1:length(cluster)
    clt.cell.idx <- split(x = cell.idx, f = cluster)
    max.probs <- apply(X = probs, MARGIN = 1, FUN = function(x) max(x)) 
    clt.probs <- split(x = max.probs, f = cluster)
    clt.batch <- split(x = batch, f = cluster)
    # Split divided idx & probs by batch
    clt.probs.batch <- lapply(X = clts, FUN = function(x) {
        split(x = clt.probs[[x]], f = clt.batch[[x]])
    })
    clt.cell.idx.batch <- lapply(X = clts, FUN = function(x) {
        split(x = clt.cell.idx[[x]], f = clt.batch[[x]])
    })
    # Quantile - median by default
    clt.b.q.split <- lapply(X = clt.probs.batch, FUN = function(x) {
        lapply(X = x, FUN = function(y) {
            quantile(y, probs = q.split) 
        })
    })
    clt.q.probs.split <- lapply(X = clts, FUN = function(x) {
        unlist(
            lapply(X = batches, FUN = function(y) {
                (clt.probs.batch[[x]][[y]] <= clt.b.q.split[[x]][[y]])
            })
        )
    })
    clt.cell.idx. <- lapply(X = clt.cell.idx.batch, FUN = function(x) unlist(x[batches])) 
    # Make new clusterings
    new.idx <- new.clt <- vector(mode="numeric", length = length(cluster))
    s.pos <- 1
    e.pos <- i <- 0
    for (clt in clts) {
        new.clt1 <- clt.cell.idx.[[clt]][clt.q.probs.split[[clt]]]
        new.clt2 <- clt.cell.idx.[[clt]][!clt.q.probs.split[[clt]]]
        make.clt <- rep(c(1+i, 2+i), c(length(new.clt1), length(new.clt2)))
        e.pos <- e.pos + length(make.clt) 
        new.idx[s.pos:e.pos] <- c(new.clt1, new.clt2)
        new.clt[s.pos:e.pos] <- make.clt
        i <- i + 2
        s.pos <- s.pos + length(make.clt) 
    }
    new.clt <- factor(new.clt[order(new.idx)])
    return(new.clt)
}