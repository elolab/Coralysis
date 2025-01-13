#' A toy dataset with 500 cells downsampled from the pbmc3k dataset.
#'
#' The preprocessing was done using Cell Ranger v2.2.0 and
#' the GRCh38.p12 human reference genome.
#' The Normalization was done using the LogNormalize method of Seurat v3 R package.
#' The sampling was done using the sample() function without replacement
#' and set.seed(1) as initialization.
#' @docType data
#'
#' @usage data(pbmc3k_500)
#'
#' @format pbmc3k_500, dgCMatrix object
#'
#' @keywords datasets
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression}
#'
#' @examples
#' data(pbmc3k_500)
"pbmc3k_500"

#' A toy SCE dataset of two 10X PBMCs 3' assays downsampled to 2K cells.
#'
#' Two 10X PBMCs 3' assays, V1 and V2, downsampled to 2K cells (1K per assay).  
#' Both datasets have been downloaded from 10X website (see source below). The 
#' PBMC dataset V1 corresponds to sample pbmc6k and V2 pbmc8k. Cells were  
#' annotated using the annotations provided by Korsunsky et al., 2019 (PMID:31740819; 
#' Source Data Figure 4 file). 
#' 
#' @docType data
#'
#' @usage 
#' data(pbmc_10Xassays)
#' pbmc_10Xassays
#'
#' @format SCE object
#'
#' @keywords datasets
#'
#' @source 
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc6k}
#' 
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k}
#'
#' @examples
#' data(pbmc_10Xassays)
#' pbmc_10Xassays
"pbmc_10Xassays"


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

#' @title Aggregates cell gene expression by clusters
#'
#' @description
#' The function aggregates cell gene expression by clusters provided.
#'
#' @param mtx Matrix with genes vs cells (rows vs cols) with gene expression to 
#' aggregate.
#' @param cluster Cluster identities vector corresponding to the cells in 
#' \code{mtx}.
#' @param select.features Should features, i.e., genes, be selected. By default 
#' \code{NULL}, all genes used. 
#' @param fun Character specifying if gene expression should be aggregated by 
#' \code{mean} or \code{sum}. By default \code{"mean"}. 
#'
#' @return Matrix of gene expressed aggregated by clusters.
#'
#' @keywords aggregated gene expression
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

#' @title Scale sparse matrix by genes (column) by batch
#'
#' @description Scales genes by batch
#'
#' @param x A matrix of class `dgCMatrix`. Cells by genes (rows x columns).  
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
#' sample's gene expression data with genes in rows and cells in columns.
#' Default is \code{NULL}.
#' @param training.ident A named factor containing sample's cluster labels for
#' each cell in training.sparse.matrix. Default is \code{NULL}.
#' @param C Cost of constraints violation in L1-regularized logistic
#' regression (C). Default is \code{0.3}.
#' @param reg.type "L1" for LASSO and "L2" for Ridge. Default is "L1".
#' @param test.sparse.matrix A sparse matrix (dgCMatrix) containing test
#' sample's gene expression data with genes in rows and cells in columns.
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

#' @title Select top or bottom N genes based on a selection criterion
#'
#' @description
#' The SelectTopGenes function enables selecting top or bottom N genes based
#' on a criterion (e.g. log2FC or adj.p.value).
#'
#' @param gene.markers A data frame of the gene markers found by
#' FindAllGeneMarkers function.
#' @param top.N How many top or bottom genes to select. Default is \code{10}.
#' @param criterion.type Which criterion to use for selecting the genes.
#' Default is "log2FC".
#' @param inverse Whether to select bottom instead of top N genes.
#' Default is \code{FALSE}.
#'
#' @return an object of `data.frame` class
#'
#' @keywords select top bottom N genes
#'
#' @importFrom dplyr group_by %>% top_n
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
#' gene_markers <- FindAllGeneMarkers(sce)
#' ## Select top 10 markers based on log2 fold-change
#' top10_log2FC <- SelectTopGenes(gene_markers,
#'                                top.N = 10,
#'                                criterion.type = "log2FC",
#'                                inverse = FALSE)
#'
#' @export
#'
SelectTopGenes <- function(gene.markers = NULL, top.N = 10,
                           criterion.type = "log2FC", inverse=FALSE)
{

  if (inverse)
  {
    top.N <- -top.N
  }

  gene.markers %>%
    group_by(.data$cluster) %>% top_n(top.N, get(criterion.type)) -> top_N

  return(top_N)
}
