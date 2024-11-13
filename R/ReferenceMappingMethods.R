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
#' @param scale.query.by Should the query data be scaled by \code{cell} or by 
#' \code{gene}. By default is \code{NULL}, i.e., is not scaled. Scale it if 
#' reference was scaled.  
#' @param project.umap Project query data onto reference UMAP (logical). By 
#' default \code{FALSE}. If \code{TRUE}, the \code{ref} object needs to have a 
#' UMAP embedding obtained with \code{RunUMAP(..., return.model = TRUE)} function. 
#' @param selecte.icp.models Select the reference ICP models to use for query 
#' cluster probability prediction. By default \code{metadata(ref)$iloreg$pca.params$select.icp.tables}, 
#' i.e., the models selected to compute the reference PCA are selected. 
#' If \code{NULL} all are used. Otherwise a numeric vector should be given
#' to select the ICP models of interest.    
#' @param k.nn The number of \code{k} nearest neighbors to use in the classification
#' KNN algorithm used to transfer labels from the reference to queries (integer).
#' By default \code{10}.  
#' @param dimred.name.prefix Dimensional reduction name prefix to add to the 
#' computed PCA and UMAP. By default nothing is added, i.e., 
#' \code{dimred.name.prefix = ""}.
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
ReferenceMapping.SingleCellExperiment <- function(ref, query, ref.label,
                                                  scale.query.by, project.umap, 
                                                  select.icp.models, k.nn, 
                                                  dimred.name.prefix) {
    # Check input params
    stopifnot(is(ref, "SingleCellExperiment"), is(query, "SingleCellExperiment"), 
              (ref.label %in% colnames(colData(ref))), any(is.null(scale.query.by), (scale.query.by %in% c("cell", "gene"))), 
              is.logical(project.umap), any(is.null(select.icp.models), is.numeric(select.icp.models)), 
              all(is.numeric(k.nn), (length(k.nn)==1)), is.character(dimred.name.prefix))
    if (is.null(metadata(ref)$iloreg$pca.model)) {
        stop("PCA model does not exist. Run 'RunPCA(...)' with 'return.model = TRUE'.")
    }
    if (project.umap & is.null(metadata(ref)$iloreg$umap.model)) {
        stop("UMAP model does not exist. Run 'RunUMAP(...)' with 'return.model = TRUE'.")
    }
    
    # Filter out genes
    ref.genes <- row.names(ref)
    query.genes <- row.names(query)
    if (!all(ref.genes %in% query.genes)) {
        message(paste0("Reference does not share all the genes with query:\n", 
                       "Reference genes: ", length(ref.genes), "\n",
                       "Query genes shared: ", sum(ref.genes %in% query.genes), 
                       "\nContinuing analysis..."))
    }
    pick.genes <- which(ref.genes %in% query.genes)
    #query <- query[ref.genes[pick.genes],]
    
    # Get model data
    if (is.null(select.icp.models)) {
        n.icps <- length(metadata(ref)$iloreg$joint.probability)
        select.icp.models <- 1:n.icps
        select.icp.tables <- metadata(ref)$iloreg$pca.params$select.icp.tables
    } else {
        select.icp.tables <- seq_along(metadata(ref)$iloreg$pca.params$select.icp.tables)
    }
    models <- metadata(ref)$iloreg$models[select.icp.models]
    pca.model <- metadata(ref)$iloreg$pca.model
    
    # Predict cluster probabilities
    query.data <- t(logcounts(query[ref.genes[pick.genes],]))
    if (!is.null(scale.query.by)) {
        if (scale.query.by=="cell") {
            query.data <- Scale(x = query.data, scale.by = "row")
        } 
        if (scale.query.by=="gene") {
            query.data <- Scale(x = query.data, scale.by = "col")
        }
    }
    colnames(query.data) <- paste0("W", pick.genes)
    query.probs <- list()
    for (m in seq_along(models)) {
        models[[m]]$W <- models[[m]]$W[, c(colnames(query.data), "Bias"), drop=FALSE] 
        pred <- predict(models[[m]], query.data, proba = TRUE)
        query.probs[[m]] <- pred$probabilities
    }
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
    reducedDim(x = query, type = paste0(dimred.name.prefix, "PCA")) <- query.pca
    query[["ilo_labels"]] <- preds.labels
    query[["ilo_labels_prob"]] <- attr(preds.labels, "prob")
    
    # Project data onto ref UMAP
    if (project.umap) {
        umap.model <- metadata(ref)$iloreg$umap.model
        if (is(umap.model, "umap")) { # from 'umap::umap' - class 'umap'
            dims <- seq_len(ncol(metadata(ref)$iloreg$umap.model$data))
            query.umap <- predict(umap.model, query.pca[,dims])
        } else { # from 'uwot::umap' - class 'list'
            dims <- seq_len(metadata(ref)$iloreg$umap.model$metric$euclidean$ndim)
            query.umap <- uwot::umap_transform(X = query.pca[,dims], model = umap.model)
        }
        row.names(query.umap) <- colnames(query)
        reducedDim(x = query, type = paste0(dimred.name.prefix, "UMAP")) <- query.umap
    }
    return(query)
}
#' @rdname ReferenceMapping
#' @aliases ReferenceMapping
setMethod("ReferenceMapping", 
          signature(ref = "SingleCellExperiment", query = "SingleCellExperiment"),
          ReferenceMapping.SingleCellExperiment)
