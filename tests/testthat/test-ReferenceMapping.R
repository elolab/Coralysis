# Packages
library("SingleCellExperiment")

set.seed(49271)
test_that("ReferenceMapping works without UMAP projection", {
  # Create data
  batches <- c("b1", "b2")
  batch <- sample(x = batches, size = nrow(iris), replace = TRUE)
  sce <- SingleCellExperiment(
    assays = list(logcounts = t(iris[, 1:4])),
    colData = DataFrame(
      "Species" = iris$Species,
      "Batch" = batch
    )
  )
  colnames(sce) <- paste0("samp", 1:ncol(sce))
  ref <- sce[, sce$Batch == "b1"]
  query <- sce[, sce$Batch == "b2"]
  
  # Run function
  ref <- RunParallelDivisiveICP(
    object = ref, k = 2, L = 25, C = 1,
    train.k.nn = 10, train.k.nn.prop = NULL,
    use.cluster.seed = FALSE,
    build.train.set = FALSE, ari.cutoff = 0.1,
    threads = 2
  )
  ref <- RunPCA(ref, p = 5, return.model = TRUE, pca.method = "stats")
  query <- ReferenceMapping(
    ref = ref, query = query, ref.label = "Species",
    project.umap = FALSE
  )
  preds_x_truth <- table(query$coral_labels, query$Species)
  
  # Test output 
  exp.true.pos <- c("setosa" = 24, "versicolor" = 24, "virginica" = 31)
  expect_equal(exp.true.pos, diag(preds_x_truth))
})

set.seed(301135)
test_that("ReferenceMapping works with UMAP projection", {
  # Create data
  batches <- c("b1", "b2")
  batch <- sample(x = batches, size = nrow(iris), replace = TRUE)
  sce <- SingleCellExperiment(
    assays = list(logcounts = t(iris[, 1:4])),
    colData = DataFrame(
      "Species" = iris$Species,
      "Batch" = batch
    )
  )
  colnames(sce) <- paste0("samp", 1:ncol(sce))
  ref <- sce[, sce$Batch == "b1"]
  query <- sce[, sce$Batch == "b2"]
  
  # Run function
  ref <- RunParallelDivisiveICP(
    object = ref, k = 2, L = 25, C = 1,
    train.k.nn = 10, train.k.nn.prop = NULL,
    use.cluster.seed = FALSE,
    build.train.set = FALSE, ari.cutoff = 0.1,
    threads = 2
  )
  ref <- RunPCA(ref, p = 5, return.model = TRUE, pca.method = "stats")
  ref <- RunUMAP(ref, return.model = TRUE)
  query <- ReferenceMapping(
    ref = ref, query = query, ref.label = "Species",
    project.umap = TRUE
  )
  preds_x_truth <- table(query$coral_labels, query$Species)
  
  # Test output 
  exp.true.pos <- c("setosa" = 31, "versicolor" = 24, "virginica" = 21)
  expect_equal(exp.true.pos, diag(preds_x_truth))
  expect_true("UMAP" %in% reducedDimNames(query))
})
