#' @title Run ICP runs parallerly with increasing number of Ks
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
#' @param k A positive integer greater or equal to \code{2}, denoting
#' the number of clusters in Iterative Clustering Projection (ICP).
#' Decreasing \code{k} leads to smaller cell populations diversity
#' and vice versa. Default is \code{15}.
#' @param d A numeric greater than \code{0} and smaller than \code{1} that
#' determines how many cells \code{n} are down- or oversampled from each cluster
#' into the training data (\code{n=N/k*d}), where \code{N} is the total number
#' of cells, \code{k} is the number of clusters in ICP. Increasing above 0.3
#' leads greadually to smaller cell populations diversity.
#' Default is \code{0.3}.
#' @param L A positive integer greater than \code{1} denoting the number of
#' the ICP runs to run. Default is \code{200}. Increasing recommended with
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
#' to randomly select for each ICP run from the complete data set. 
#' This is a new feature intended to speed up the process
#' with larger data sets. Default is \code{Inf}, which means using all cells.
#' @param train.with.bnn Train data with batch nearest neighbors. Default is 
#' \code{TRUE}. Only used if \code{batch.label} is given.   
#' @param train.k.nn Train data with batch nearest neighbors using \code{k} 
#' nearest neighbors. Default is \code{10}. Only used if \code{train.with.bnn} 
#' is \code{TRUE}.
#' @param build.train.set Logical specifying if a training set should be built 
#' from the data or the whole data should be used for training. By default 
#' \code{FALSE}.
#' @param build.train.params A list of parameters to be passed to the function
#' \code{AggregateDataByBatch()}.
#' @param scale A logical specifying if data should be scaled before training. 
#' Default is \code{FALSE}.
#' @param verbose A logical value to print verbose during the ICP run in case 
#' of parallelization, i.e., 'threads' different than \code{1}. Default 'FALSE'. 
#'
#' @name IntegrateData
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
IntegrateData.SingleCellExperiment <- function(object, batch.label, 
                                               k, d, L, r, C,
                                               reg.type, max.iter,
                                               threads, icp.batch.size, 
                                               train.with.bnn, train.k.nn,
                                               build.train.set, 
                                               build.train.params,
                                               scale,
                                               verbose) {
    
    if (!is(object,"SingleCellExperiment")) {
        stop("object must of 'sce' class")
        return(object)
    }
    
    if (!is.numeric(k) | all(k < 2) | all(k%%1 != 0))
    {
        stop("k must be positive integer(s) and greater than 1")
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
    
    if (build.train.set) {
        build.train.params <- c(list(object = object, batch.label = batch.label), build.train.params)
        clustered.object <- do.call(AggregateDataByBatch, build.train.params)
        dataset <- t(logcounts(clustered.object))
    } else {
        dataset <- t(logcounts(object))
    }
    
    if (!is.null(batch.label)) {
        if (build.train.set) {
            batch.label <- as.character(clustered.object[["batch"]])
            names(batch.label) <- colnames(clustered.object)
        } else {
            batch.label <- as.character(object[[batch.label]])
            names(batch.label) <- colnames(object)
        }
    }
    
    if (scale) {
        dataset <- Scale(x = as(dataset, "sparseMatrix"), scale.by="row")
    }
    parallelism <- TRUE
    
    icp <- list()
    for (clt in k) {
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
        k.clt <- paste0("k", clt)
        if (parallelism) {
            pb <- txtProgressBar(min = 1, max = L, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
            icp[[k.clt]] <- foreach(task = seq_len(L),
                                    .verbose = FALSE,
                                    .combine = list,
                                    .maxcombine = 1000,
                                    .inorder = FALSE,
                                    .multicombine = TRUE,
                                    .packages=c("ILoReg2", "parallel"),
                                    .options.snow = opts)  %dorng% {
                                        tryCatch(expr = {
                                            message(paste0("\nICP run: ", task))
                                            RunICP(normalized.data = dataset, batch.label = batch.label, 
                                                   k = clt, d = d, r = r, C = C, reg.type = reg.type, 
                                                   max.iter = max.iter, icp.batch.size = icp.batch.size, 
                                                   train.with.bnn = train.with.bnn, train.k.nn = train.k.nn)
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
            for (l in seq_len(L)) {
                try({
                    message(paste0("ICP run: ",l))
                    res <- RunICP(normalized.data = dataset, batch.label = batch.label, 
                                  k = clt, d = d, r = r, C = C, reg.type = reg.type, 
                                  max.iter = max.iter, icp.batch.size = icp.batch.size, 
                                  train.with.bnn = train.with.bnn, train.k.nn = train.k.nn)
                    icp[[k.clt]][[l]] <- res
                })
            }
        }
    }
    out <- do.call(c, icp)
    
    metadata(object)$iloreg$metrics <-
        lapply(out, function(x) x$metrics)
    
    metadata(object)$iloreg$models <-
        lapply(out, function(x) x$model)
    
    if (build.train.set) {
        test.data <- t(logcounts(object))
        colnames(test.data) <- paste0("W", 1:ncol(test.data))
        if (scale) {
            test.data <- Scale(x = as(test.data, "sparseMatrix"), scale.by="row")
        }
        metadata(object)$iloreg$joint.probability <- 
            lapply(metadata(object)$iloreg$models, function(x) {
                predict(x, test.data, proba=TRUE)$probabilities
            })
    } else {
        metadata(object)$iloreg$joint.probability <-
            lapply(out, function(x) x$probabilities)
    }
    
    # Order output lists by increasing standard deviation of cluster probability tables 
    icp.k <- rep(k, each = L)
    icp.k.idx <- split(x = 1:length(icp.k), f = icp.k)
    order.list <- lapply(X = icp.k.idx, FUN = function(x) {
        order(unlist(lapply(metadata(object)$iloreg$joint.probability[x], sd))) 
        })
    ordered.list <- unlist(lapply(X = 1:length(k), FUN = function(x) {
        icp.k.idx[[x]][order.list[[x]]]
    }))
    metadata(object)$iloreg$joint.probability <- metadata(object)$iloreg$joint.probability[ordered.list]
    metadata(object)$iloreg$metrics <- metadata(object)$iloreg$metrics[ordered.list]
    metadata(object)$iloreg$models <- metadata(object)$iloreg$models[ordered.list]
    
    return(object)
}
#' @rdname IntegrateData
#' @aliases IntegrateData
setMethod("IntegrateData", signature(object = "SingleCellExperiment"),
          IntegrateData.SingleCellExperiment)
