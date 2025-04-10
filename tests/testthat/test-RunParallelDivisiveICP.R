# Packages
library("SingleCellExperiment")

set.seed(951355)
test_that("RunParallelDivisiveICP works with parallelization", {
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
  
  # Run function
  sce <- RunParallelDivisiveICP(
    object = sce, batch.label = "Batch",
    k = 2, L = 25, C = 1, train.k.nn = 10,
    train.k.nn.prop = NULL, use.cluster.seed = FALSE,
    build.train.set = FALSE, ari.cutoff = 0.1,
    threads = 2
  )
  sce <- SummariseCellClusterProbability(
    object = sce, icp.run = 13, 
    icp.round = 1, save.in.sce = TRUE
  )
  
  # Test output
  exp.clt <- c(
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", 
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2",
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2",
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2"
  )
  expect_equal(25*log2(2), length(GetCellClusterProbability(sce, concatenate = FALSE)))
  expect_equal(exp.clt, sce$icp_run_round_13_1_clusters)
})

set.seed(124056)
test_that("RunParallelDivisiveICP works without parallelization", {
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
  
  # Run function
  sce <- RunParallelDivisiveICP(
    object = sce, batch.label = "Batch",
    k = 2, L = 25, C = 1, train.k.nn = 10,
    train.k.nn.prop = NULL, use.cluster.seed = FALSE,
    build.train.set = FALSE, ari.cutoff = 0.1,
    threads = 1
  )
  sce <- SummariseCellClusterProbability(
    object = sce, icp.run = 5, 
    icp.round = 1, save.in.sce = TRUE
  )
  
  # Test output
  exp.clt <- c(
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", 
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", 
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2"
  )
  expect_equal(25*log2(2), length(GetCellClusterProbability(sce, concatenate = FALSE)))
  expect_equal(exp.clt, sce$icp_run_round_5_1_clusters)
})

set.seed(126599)
test_that("RunParallelDivisiveICP works k>2", {
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
  
  # Run function
  sce <- RunParallelDivisiveICP(
    object = sce, batch.label = "Batch",
    k = 4, L = 10, C = 1, train.k.nn = 5,
    train.k.nn.prop = NULL, use.cluster.seed = FALSE,
    build.train.set = FALSE, ari.cutoff = 0.1,
    threads = 2
  )
  sce <- SummariseCellClusterProbability(
    object = sce, icp.run = NULL, 
    icp.round = NULL, save.in.sce = TRUE
  )
  
  # Test output
  exp.clt <- c(
    "4", "3", "4", "3", "4", "4", "4", "4", "3", "3", "4", "3", "3", "4", "4", "4", "4", "4",
    "4", "4", "3", "4", "4", "4", "3", "3", "4", "4", "4", "3", "3", "4", "4", "4", "3", "4", 
    "4", "4", "4", "4", "4", "3", "4", "4", "4", "4", "4", "4", "4", "4", "1", "1", "2", "2", 
    "2", "2", "1", "1", "2", "1", "2", "1", "2", "2", "1", "1", "1", "1", "2", "2", "1", "1", 
    "2", "2", "1", "1", "2", "2", "2", "1", "2", "1", "1", "2", "1", "1", "1", "2", "1", "2", 
    "2", "2", "1", "1", "2", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2"
  )
  expect_equal(10*log2(4), length(GetCellClusterProbability(sce, concatenate = FALSE)))
  expect_equal(exp.clt, sce$icp_run_round_10_2_clusters)
})

set.seed(86661)
test_that("RunParallelDivisiveICP works without a batch label", {
  # Create data
  sce <- SingleCellExperiment(
    assays = list(logcounts = t(iris[, 1:4])),
    colData = DataFrame(
      "Species" = iris$Species
    )
  )
  colnames(sce) <- paste0("samp", 1:ncol(sce))
  
  # Run function
  sce <- RunParallelDivisiveICP(
    object = sce, batch.label = NULL,
    k = 2, L = 25, C = 1, train.k.nn = 10,
    train.k.nn.prop = NULL, use.cluster.seed = FALSE,
    build.train.set = FALSE, ari.cutoff = 0.1,
    threads = 2
  )
  sce <- SummariseCellClusterProbability(
    object = sce, icp.run = 21, 
    icp.round = 1, save.in.sce = TRUE
  )
  
  # Test output
  exp.clt <- c(
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
    "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
    "2", "2", "2", "2", "2", "2"
  )
  expect_equal(25*log2(2), length(GetCellClusterProbability(sce, concatenate = FALSE)))
  expect_equal(exp.clt, sce$icp_run_round_21_1_clusters)
})
