#' @title Reference mapping
#'
#' @description
#' This functions allows to project new query data sets onto a reference built 
#' with ILoReg. 
#'
#' @param ref An object of \code{SingleCellExperiment} class trained with ILoReg.
#' @param query An object of \code{SingleCellExperiment} class to project onto 
#' \code{ref}.
#' @param scale.query Should the query logcounts be scaled or not (logical). By 
#' default \code{TRUE}. Scale it if reference was scaled.  
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
#' 
ReferenceMapping.SingleCellExperiment <- function(ref, query, scale.query) {
    # Filter out genes
    ref.genes <- row.names(ref)
    query.genes <- row.names(query)
    pick.genes <- which(ref.genes %in% query.genes)
    query <- query[ref.genes[pick.genes],]
    
    # Get model data
    models <- metadata(ref)$iloreg$models
    pca.model <- metadata(ref)$iloreg$pca.model
    lda.model <- metadata(ref)$iloreg$lda.model
    
    # Predict
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
    probs <- list("joint.probability" = list("ref" = metadata(ref)$iloreg$joint.probability, 
                                             "query" = query.probs))
    query.probs <- do.call(cbind, query.probs)
    query.pca <- scale(query.probs, pca.model$center, pca.model$scale) %*% pca.model$rotation 
    query.lda <- predict(lda.model, query.pca)
    
    # Object
    ref.coldata <- colData(ref)
    query.coldata <- colData(query)
    shared.cols <- intersect(colnames(ref.coldata), colnames(query.coldata))
    ref.uniq.cols <- colnames(ref.coldata)[!colnames(ref.coldata) %in% shared.cols]  
    query.uniq.cols <- colnames(query.coldata)[!colnames(query.coldata) %in% shared.cols]
    all.cols <- c(shared.cols, ref.uniq.cols, query.uniq.cols)
    ref.coldata[,query.uniq.cols] <- NA
    query.coldata[,ref.uniq.cols] <- NA
    coldata <- rbind(ref.coldata[,all.cols], query.coldata[,all.cols])
    pca <- rbind(reducedDim(ref, "PCA"), query.pca)
    lda <- rbind(reducedDim(ref, "LDA"), query.lda$x)
    row.names(pca) <- row.names(lda) <- row.names(coldata)
    
    log.counts <- cbind(logcounts(ref[pick.genes,]), logcounts(query))
    ref.map <- SingleCellExperiment(
        assays = list("logcounts" = log.counts), 
        colData = coldata,
        reducedDims = list("PCA" = pca, "LDA" = lda),
        metadata = list("iloreg" = probs)
    )
    return(ref.map)
}
#' @rdname ReferenceMapping
#' @aliases ReferenceMapping
setMethod("ReferenceMapping", 
          signature(ref = "SingleCellExperiment", query = "SingleCellExperiment"),
          ReferenceMapping.SingleCellExperiment)
