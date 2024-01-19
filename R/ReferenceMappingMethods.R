#' @title Reference mapping
#'
#' @description
#' This function allows to project new query data sets onto a reference built 
#' with ILoReg as well as transfer cell labels from the reference to queries. 
#'
#' @param ref An object of \code{SingleCellExperiment} class trained with ILoReg
#' and after running \code{RunPCA(..., return.model = TRUE)} function. 
#' @param query An object of \code{SingleCellExperiment} class to project onto 
#' \code{ref}.
#' @param ref.label A character cell metadata column name from the \code{ref} 
#' object to transfer to the queries.  
#' @param scale.query Should the query logcounts be scaled or not (logical). By 
#' default \code{TRUE}. Scale it if reference was scaled.  
#' @param project.umap Project query data onto reference UMAP (logical). By 
#' default \code{FALSE}. If \code{TRUE}, the \code{ref} object needs to have a 
#' UMAP embedding obtained with \code{RunUMAP(..., return.model = TRUE)} function. 
#' @param selecte.icp.models Select the reference ICP models to use for query 
#' cluster probability prediction. By default \code{NULL}, i.e., all are used.
#' A vector of \code{integers} should be given otherwise.   
#' @param k.nn The number of \code{k} nearest neighbors to use in the classification
#' KNN algorithm used to transfer labels from the reference to queries (integer).
#' By default \code{10}.  
#'
#' @name ReferenceMapping
#'
#' @return An object of \code{SingleCellExperiment} class.
#'
#' @keywords iterative clustering projection ICP logistic regression LIBLINEAR
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment logcounts SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom class knn
#' 
ReferenceMapping.SingleCellExperiment <- function(ref, query, ref.label,
                                                  scale.query, project.umap, 
                                                  select.icp.models, k.nn) {
    # Check input params
    if (is.null(metadata(ref)$iloreg$pca.model)) {
        stop("PCA model does not exist. Run 'RunPCA(...)' with 'return.model = TRUE'.")
    }
    if (project.umap & is.null(metadata(ref)$iloreg$umap.model)) {
        stop("UMAP model does not exist. Run 'RunUMAP(...)' with 'return.model = TRUE'.")
    }

    # Filter out genes
    ref.genes <- row.names(ref)
    query.genes <- row.names(query)
    pick.genes <- which(ref.genes %in% query.genes)
    query <- query[ref.genes[pick.genes],]
    
    # Get model data
    if (is.null(select.icp.models)) {
        n.icps <- length(metadata(ref)$iloreg$joint.probability)
        select.icp.models <- 1:n.icps
    }
    models <- metadata(ref)$iloreg$models[select.icp.models]
    pca.model <- metadata(ref)$iloreg$pca.model
    
    # Predict cluster probabilities
    query.data <- t(logcounts(query))
    if (scale.query) {
        query.data <- Scale(x = query.data, scale.by = "row")
    }
    colnames(query.data) <- paste0("W", pick.genes)
    query.probs <- list()
    for (m in seq_along(models)) {
        models[[m]]$W <- models[[m]]$W[, c(colnames(query.data), "Bias"), drop=FALSE] 
        pred <- predict(models[[m]], query.data, proba = TRUE)
        query.probs[[m]] <- pred$probabilities
    }
    select.icp.tables <- metadata(ref)$iloreg$pca.params$select.icp.tables
    probs <- do.call(cbind, query.probs[select.icp.tables])
    
    # Project data onto ref PCA
    query.pca <- scale(probs, pca.model$center, pca.model$scale) %*% pca.model$rotation 
    
    # Predict classifications with KNN
    preds.labels <- class::knn(train = pca.model$x, test = query.pca, 
                               cl = ref[[ref.label]], k = k.nn, 
                               prob = TRUE)
    
    # Add predictions to query object
    metadata(query)$iloreg <- list()
    metadata(query)$iloreg$joint.probability <- query.probs   
    reducedDim(query, type = "PCA") <- query.pca
    query[["ilo_labels"]] <- preds.labels
    query[["ilo_labels_prob"]] <- attr(preds.labels, "prob")
    
    # Project data onto ref UMAP
    if (project.umap) {
        umap.model <- metadata(ref)$iloreg$umap.model
        query.umap <- predict(umap.model, query.pca)
        row.names(query.umap) <- colnames(query)
        reducedDim(query, "UMAP") <- query.umap
    }
    return(query)
}
#' @rdname ReferenceMapping
#' @aliases ReferenceMapping
setMethod("ReferenceMapping", 
          signature(ref = "SingleCellExperiment", query = "SingleCellExperiment"),
          ReferenceMapping.SingleCellExperiment)
