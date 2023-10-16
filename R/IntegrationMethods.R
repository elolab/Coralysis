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
#' @param k A positive integer greater or equal to \code{2}, denoting
#' the number of clusters in Iterative Clustering Projection (ICP).
#' Decreasing \code{k} leads to smaller cell populations diversity
#' and vice versa. Default is \code{8}.
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
                res <- RunDivisiveICP(normalized.data = dataset, batch.label = batch.label, 
                                      k = k, d = d, r = r, C = C, reg.type = reg.type, 
                                      max.iter = max.iter, icp.batch.size = icp.batch.size, 
                                      train.with.bnn = train.with.bnn, train.k.nn = train.k.nn)
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
    metadata(object)$iloreg$joint.probability <- unlist(metadata(object)$iloreg$joint.probability, recursive = FALSE)
    
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
#' @param k A positive integer greater or equal to 2, denoting the number of
#' clusters in ICP. Default is \code{8}.
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
                           k = 8, d = 0.3, r = 5, C = 5,
                           reg.type = "L1", max.iter = 200, 
                           icp.batch.size=Inf, train.with.bnn = TRUE, 
                           train.k.nn = 10) {
    
    #first_round <- TRUE
    metrics <- NULL
    idents <- list()
    iterations <- 1
    #probs <- NULL
    
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
    for (k in Ks) {
        first_round <- TRUE
        i <- i + 1
        while (TRUE) {
            
            # Step 1: initialize clustering (ident_1) randomly, ARI=0 and r=0
            if (first_round) {
                if (k == 2) {
                    ident_1 <- factor(sample(seq_len(k), nrow(normalized.data), replace = TRUE))
                } else {
                    ident_1 <- RandomlyDivisiveClustering(cluster=preds[[i-1]], k=2)
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
                                        k = train.k.nn)
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
            if (length(levels(factor(as.character(ident_2)))) < k) {
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
                if (ari>0.5) {
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
