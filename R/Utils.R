#' @title Cluster cells 
#'
#' @description
#' The function clusters cells with the K-means++ algorithm
#'
#' @param object An object of \code{SingleCellExperiment} class.
#' @param nclusters Cluster the cells into n clusters. Ignored if the number of 
#' cells in \code{object} is lower or equal to \code{nclusters}. 
#' @param use.emb Should the embedding be used to cluster or the log-transformed 
#' data. By default \code{TRUE}. 
#' @param emb.name Which embedding to use. By default \code{"PCA"}. 
#'
#' @return A SingleCellExperiment object with clusters.
#'
#' @keywords cluster
#' 
#' @importFrom flexclust kcca kccaFamily
#'
ClusterCells <- function(object, nclusters=500, use.emb=TRUE, emb.name="PCA") {
    
    # Check input params
    stopifnot(is(object, "SingleCellExperiment"), 
              (is.numeric(nclusters) && (nclusters%%1 == 0)), 
              (is.logical(use.emb)))
    
    if (use.emb) {
        stopifnot(is.character(emb.name) && (length(emb.name)==1) && (emb.name %in% reducedDimNames(object)))
        data.cluster <- reducedDim(object, emb.name)
    } else {
        data.cluster <- t(logcounts(object))
    }
    no.cells <- nrow(data.cluster)
    if (no.cells <= nclusters) {
        message("No. of cells is lower or equal than the no. of 'nclusters': ", 
                no.cells, " <= ", nclusters, "\nAll the cells will be included and the clustering step will be skipped!\n")
        nclusters <- no.cells
    }
    clusters <- kcca(x = data.cluster, k = nclusters, family = kccaFamily("kmeans"), 
                     control = list(initcent = "kmeanspp"))
    metadata(object)$clusters <- clusters 

    return(object)
} 

#' @title Aggregates cell feature expression by clusters
#'
#' @description The function aggregates cell feature expression by clusters provided.
#'
#' @param mtx Matrix with features vs cells (rows vs cols) with feature expression to 
#' aggregate.
#' @param cluster Cluster identities vector corresponding to the cells in 
#' \code{mtx}.
#' @param select.features Should features be selected. By default \code{NULL}, 
#' all features used. 
#' @param fun Character specifying if feature expression should be aggregated by 
#' \code{mean} or \code{sum}. By default \code{"mean"}. 
#'
#' @return Matrix of feature expressed aggregated by clusters.
#'
#' @keywords aggregated feature expression
#'
AggregateClusterExpression <- function(mtx, cluster, select.features = NULL, 
                                       fun = "mean") {
    stopifnot(fun %in% c("sum", "mean"))
    # Convert mtx to df
    if(is(mtx, "dgCMatrix")) mtx <- as.data.frame(as.matrix(mtx))
    if (is.matrix(mtx)) mtx <- as.data.frame(mtx)
    if (!is.null(select.features)) mtx <- mtx[select.features,]
    # Calculate average expression by cluster
    #based on: https://github.com/neurorestore/Libra/blob/main/R/to_pseudobulk.R#L87
    mm <- model.matrix(~ 0 + cluster, data=mtx)
    ps <- as.matrix(mtx) %*% mm
    if (fun=="mean") {
        scale_fct <- as.vector(table(cluster))
        ps <- sweep(ps, 2, scale_fct, "/")
    }
    return(ps)
}


#' @title Scale a sparse matrix by row or column
#'
#' @description Faster implementation of \code{scale} function. It diverges from 
#' the \code{scale} function by not performing the root-mean-square when 
#' \code{scale=TRUE} and \code{center=FALSE}. In this case it divides the values
#' by the standard deviation of the column or row used (depending on \code{scale.by}). 
#'
#' @param x A matrix of class `dgCMatrix`. 
#' @param center A logical. By default \code{TRUE}. Subtract the values by the 
#' row or column mean depending on the `scale.by` parameter. 
#' @param scale A logical. By default \code{TRUE}. Divide the values by the row 
#' or column standard deviation depending on the `scale.by` parameter
#' @param scale.by Scale by `row` or `col` (=column), i.e., use the row or 
#' column mean and/or standard deviations to center and /or scale the data. 
#' Default is \code{col}.
#' 
#' @return A matrix of class `dgCMatrix`. 
#'
#' @keywords scale
#'
#' @importFrom sparseMatrixStats colMeans2 colSds rowMeans2 rowSds
#' 
Scale <- function(x, center=TRUE, scale=TRUE, scale.by="col") {
    # Check input
    stopifnot(is.logical(center) & is.logical(scale))
    stopifnot(inherits(x, what="dgCMatrix"))
    stopifnot(scale.by %in% c("col", "row"))
    # Assign funs 
    col.mean <- sparseMatrixStats::colMeans2
    col.sds <- sparseMatrixStats::colSds
    row.mean <- sparseMatrixStats::rowMeans2
    row.sds <- sparseMatrixStats::rowSds
    # Scale
    if (center & scale & (scale.by=="col")) {
        m <- col.mean(x)
        s <- col.sds(x)
        names(m) <- names(s) <- colnames(x)
        x <- t((t(x) - m) / s) 
    } else if (center & scale & (scale.by=="row")) {
        m <- row.mean(x)  
        s <- row.sds(x)
        names(m) <- names(s) <- row.names(x)
        x <- ((x - m) / s)
    } else if (center & (scale.by=="col")) {
        m <- col.mean(x)
        names(m) <- colnames(x)
        x <- t(t(x) - m)
    } else if (center & (scale.by=="row")) {
        m <- row.mean(x)
        names(m) <- row.names(x)
        x <- x - m
    } else if (scale & (scale.by=="col")) {
        s <- col.sds(x)
        names(s) <- colnames(x)
        x <- t(t(x) / s)
    } else {
        s <- row.sds(x)
        names(s) <- row.names(x)
        x <- x / s
    } 
    # Convert to sparse format
    x <- Matrix::Matrix(x, sparse=T)
    return(x)
}


#' @title Scale sparse matrix by features (column) by batch
#'
#' @description Scales features by batch
#'
#' @param x A matrix of class `dgCMatrix`. Cells by features (rows x columns).  
#' @param batch A character vector with batch labels corresponding to the cells
#' given in \code{x}. The character batch labels need to be named
#' with the cells names given in the rows of \code{x}. 
#' 
#' @return A scaled matrix of class `dgCMatrix`. 
#'
#' @keywords scale
#' 
ScaleByBatch <- function(x, batch) {
    x.batch <- split(x = as.data.frame(as.matrix(x)), f = batch)
    x.batch <- lapply(X = x.batch, function(x) {
        Scale(x = as(as.matrix(x), "sparseMatrix"), scale.by = "col")
    })
    x <- do.call(rbind, x.batch)
    x <- x[names(batch),]
    x[is.na(x)] <- 0
    return(x)
}


#' @title Down- and oversample data
#'
#' @description
#' The function implements a script down- and oversamples data to
#' include n cells.
#'
#' @param x A character or numeric vector of data to down-and oversample.
#' @param n How many cells to include per cluster.
#'
#' @return a list containing the output of the LiblineaR prediction
#'
#' @keywords downsampling oversampling
#'
DownOverSampling <- function(x, n = 50) {
  if (length(x) < n) {
    res <- sample(x, size = n, replace = TRUE)
  } else {
    res <- sample(x, size = n, replace = FALSE)
  }
  return(res)
}


#' @title Down- and oversample data evenly batches
#'
#' @description
#' The function down- and over-samples cluster cells evenly by batch.
#'
#' @param x A character or numeric vector of data to down-and oversample.
#' @param batch A character vector with batch labels corresponding to \code{x}.
#' @param n How many cells to include per cluster.
#'
#' @return a list containing the output of the LiblineaR prediction
#'
#' @keywords downsampling oversampling
#'
DownOverSampleEvenlyBatches <- function(x, batch, n = 50) {
    n.batch <- length(unique(as.character(batch)))
    sample.per.batch <- ceiling(n / n.batch)
    group.batch <- split(x = x, f = batch)
    downsample.batch <- lapply(X = group.batch, FUN = function(x) DownOverSampling(x = x, n = sample.per.batch))
    downsample.batch <- unlist(downsample.batch)
    return(downsample.batch)
}


#' @title Find batch k nearest neighbors  
#'
#' @description
#' The function finds batch k nearest neighbors for the cell with the highest 
#' probability for every batch.
#'
#' @param idx A numeric vector with the cell ids to retrieve from the data set.
#' @param group A character vector with batch labels corresponding to \code{idx}.
#' @param prob A numeric vector with cell probabilities corresponding to \code{idx}.
#' @param k The number of nearest neighbors to search for. Default is \code{10}.
#'
#' @return a list containing the k nearest neighbors for every cell queried
#'
#' @keywords knn 
#'
#' @importFrom RANN nn2
#'
FindBatchKNN <- function(idx, group, prob, k = 10) {
    prob.groups <- split(prob, group)
    prob.idx <- split(idx, group)
    uniq.groups <- names(prob.groups)
    comb.groups <- expand.grid(uniq.groups, uniq.groups)
    comb.groups <- comb.groups[comb.groups$Var1 != comb.groups$Var2,]
    comb.groups <- data.frame(t(comb.groups))
    obs.groups <- unlist(lapply(X = prob.idx, FUN = function(x) length(x)))
    knn.idx <- list()
    j <- 0
    for (pair in comb.groups) {
        j <- j + 1
        if (length(pair) > 0) {
            if (obs.groups[pair[1]] >= k) {
                knn <- nn2(data = prob.groups[[pair[1]]], query = prob.groups[[pair[2]]][which.max(prob.groups[[pair[2]]])], 
                           k = k)
                tmp.knn.idx <- prob.idx[[comb.groups[1, j]]][knn$nn.idx]
            } else {
                tmp.knn.idx <- prob.idx[[comb.groups[1, j]]]
            }
        } else {
            if (length(idx) >= k) {
                tmp.knn.idx <- sample(x = idx, size = k)
            } else {
                tmp.knn.idx <- sample(x = idx, size = k, replace = TRUE)
            }
        }
        knn.idx[[paste(comb.groups[, j], collapse = "_")]] <- tmp.knn.idx
    }
    return(knn.idx)
}


#' @title Find batch k nearest neighbors per cluster 
#'
#' @description
#' The function finds batch k nearest neighbors for the cell with the highest 
#' probability for every batch per cluster.
#'
#' @param preds A numeric vector with the cell cluster predictions.
#' @param probs A numeric matrix with cell cluster probabilities.
#' @param batch A character with batch labels.
#' @param k The number of nearest neighbors to search for. Default is \code{10}.
#' @param k.prop A numeric (higher than 0 and lower than 1) corresponding to the 
#' fraction of cells per cluster to use as \code{k} nearest neighbors. Default 
#' is \code{NULL} meaning that the number of \code{k} nearest neighbors is equal
#' to \code{k}. If given, \code{k} parameter is ignored and \code{k} is calculated
#' based on \code{k.prop}.  
#'
#' @return a list containing the k nearest neighbors for every cluster
#'
#' @keywords knn 
#'
FindClusterBatchKNN <- function(preds, probs, batch, k = 10, k.prop = NULL) {
    clts <- ncol(probs)
    clt.knn <- list()
    nbatches <- length(unique(as.character(batch)))
    for (clt in 1:clts) {
        idx <- which(preds == clt)
        group <- batch[idx]
        prob <- probs[idx, clt]
        if (!is.null(k.prop)) {
            k <- ceiling(length(idx) * k.prop / nbatches)
        }
        tmp <- FindBatchKNN(idx = idx, group = group, prob = prob, k = k)
        clt.knn[[clt]] <- unique(unlist(tmp))
    }
    return(clt.knn)
}


#' @title Clustering projection using logistic regression from
#' the LiblineaR R package
#'
#' @description
#' The function implements a script that downsamples data a dataset, trains
#' a logistic regression classifier model
#' and then projects its clustering onto itself using a trained
#' L1-regularized logistic regression model.
#'
#' @param training.sparse.matrix A sparse matrix (dgCMatrix) containing training
#' sample's feature expression data with features in rows and cells in columns.
#' Default is \code{NULL}.
#' @param training.ident A named factor containing sample's cluster labels for
#' each cell in training.sparse.matrix. Default is \code{NULL}.
#' @param C Cost of constraints violation in L1-regularized logistic
#' regression (C). Default is \code{0.3}.
#' @param reg.type "L1" for LASSO and "L2" for Ridge. Default is "L1".
#' @param test.sparse.matrix A sparse matrix (dgCMatrix) containing test
#' sample's feature expression data with features in rows and cells in columns.
#' Default is \code{NULL}.
#' @param d A numeric smaller than \code{1} and greater than \code{0}
#' that determines how many cells per cluster should be
#' down- and oversampled (d in N/k*d), where N is the total number of cells
#' and k the number of clusters. Default is \code{0.3}.
#' @param batch.label A character vector with batch labels corresponding to the cells
#' given in \code{training.ident}. The character batch labels need to be named
#' with the cells names given in \code{training.ident}. By default \code{NULL}, 
#' i.e., cells are sampled evenly regardless their batch. 
#' @param training_ident_subset A character or numeric vector with cell ids to 
#' use as train set. By default \code{NULL}. If given, the down- and oversampled
#' parameters are ignored.
#'
#' @return a list containing the output of the LiblineaR prediction
#'
#' @keywords logistic regression LiblineaR projection downsampling oversampling
#'
#' @import Matrix
#' @import SparseM
#' @importFrom methods as
#' @importFrom LiblineaR LiblineaR
#' @importFrom stats predict
#'
LogisticRegression <- function(training.sparse.matrix = NULL,
                               training.ident = NULL,
                               C = 0.3,
                               reg.type = "L1",
                               test.sparse.matrix = NULL,
                               d = 0.3, 
                               batch.label = NULL, 
                               training_ident_subset = NULL) {
    
    # Downsample training data
    if (!is.null(d)) {
        cells_per_cluster <- ceiling((length(training.ident) / (length(levels(training.ident)))) * d)
        if (is.null(training_ident_subset)) {
            if (is.null(batch.label)) {
                training_ident_subset <- as.character(unlist(lapply(split(names(training.ident),training.ident), function(x) DownOverSampling(x, cells_per_cluster))))
            } else {
                training_ident_subset <- as.character(unlist(lapply(split(names(training.ident), training.ident), function(x) {
                    DownOverSampleEvenlyBatches(x, batch = batch.label[x], cells_per_cluster)
                })))
            } 
        } 
        training.ident <- training.ident[training_ident_subset]
        training.sparse.matrix <- training.sparse.matrix[training_ident_subset,]
    }

  # Transform training and test data from dgCMatrix to matrix.csr
  training.sparse.matrix <- as(training.sparse.matrix, "matrix.csr")
  test.sparse.matrix <- as(test.sparse.matrix, "matrix.csr")

  if (reg.type=="L2")
  {
    type <- 7
  } else if (reg.type=="L1")
  {
    type <- 6 #L1
  } else {
    stop("'reg.type' must be either 'L1' or 'L2'")
  }

  model <- LiblineaR(training.sparse.matrix, training.ident, type = type, cost = C)
  prediction <- predict(model, proba = TRUE, test.sparse.matrix)

  return(list(prediction = prediction, model = model))
}


#' @title Divisive Iterative Clustering Projection (ICP) clustering
#'
#' @description
#' The function implements divisive Iterative Clustering Projection (ICP) clustering: 
#' a supervised learning-based clustering, which maximizes clustering similarity
#' between the clustering and its projection by logistic regression, doing it in 
#' a divisive clustering manner.
#'
#' @param normalized.data A sparse matrix (dgCMatrix) containing
#' normalized feature expression data with cells in rows and features in columns.
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
            
            message(paste0("probability matrix dimensions = ",paste(dim(res_prediction$probabilities), collapse = " ")))
            
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
#' should represent cells and features, respectively.   
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
#' @param center Should the features given in \code{data} centered before
#' performing the PCA (logical). By default \code{TRUE}. 
#' @param scale. Should the features given in \code{data} scaled before
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
    features.sd <- colSds(data) 
    features.sd.diff0 <- which(features.sd != 0)
    cell.names <- row.names(data)
    pca <- prcomp_irlba(x=data[,features.sd.diff0], n=p, center = center, scale. = scale.)
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
