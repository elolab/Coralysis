seurat_object <- readRDS("~/proj1_bm_ILoReg_int_test.rds")
    
library(Seurat)

so.list <- SplitObject(seurat_object, split.by = "Donor")

# normalize and identify variable features for each dataset independently
so.list <- lapply(X = so.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = so.list)


### scale data????
# so.list <- lapply(X = so.list, FUN = function(x) {
#     x <- ScaleData(x, features = features, verbose = FALSE)
# })

X <- seurat_object@assays$RNA@data[features,]
# X <- seurat_object@assays$RNA@data
ids <- seurat_object@meta.data$Donor

detach("package:Seurat")

suppressMessages(library(ILoReg2))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(cowplot))
# The dataset was normalized using the LogNormalize method from the Seurat R package.
set.seed(1)

probability_list <- list()
for (id in names(table(ids)))
{
    X_id <- X[,ids==id]
    
    sce <- SingleCellExperiment(assays = list(logcounts = X_id))
    sce <- PrepareILoReg2(sce)
    
    sce <- RunParallelICP(object = sce, k = 15,
                          d = 0.5, L = 100, 
                          r = 5, C = 1.0, 
                          reg.type = "L1", threads = 2)
    
    X_ <- t(X[rownames(sce),])
    colnames(X_) <- paste0("W",1:ncol(X_))
    
    # L = 100
    for (j in 1:100)
    {
        p<-predict(sce@metadata$iloreg$models[[j]],X_,proba = TRUE)
        probability_list[[paste0(id,"_",j)]] <- p$probabilities 
    }
    
}

saveRDS(probability_list,file="~/proj1_bm_ILoReg_int_test_icp_probs.rds")

sce <- SingleCellExperiment(assays = list(logcounts = X))
sce <- PrepareILoReg2(sce)

sce@metadata$iloreg$joint.probability <- probability_list

X <- scale(do.call(cbind,probability_list),center = TRUE,scale = TRUE)

A <- RSpectra::svds(X,k=50)

reducedDim(sce,"PCA") <- A$u

sce <- RunUMAP(sce)

AnnotationScatterPlot(sce,annotation = seurat_object@meta.data$Purified_cell_type,dim.reduction.type = "umap")
AnnotationScatterPlot(sce,annotation = seurat_object@meta.data$annotated_cell_identity.ontology_label,dim.reduction.type = "umap")
AnnotationScatterPlot(sce,annotation = seurat_object@meta.data$Donor,dim.reduction.type = "umap")
AnnotationScatterPlot(sce,annotation = seurat_object@meta.data$annotated_cell_identity.text,dim.reduction.type = "umap")
AnnotationScatterPlot(sce,annotation = seurat_object@meta.data$Purified_cell_type,dim.reduction.type = "umap")


load("C:/Users/pajosm/AppData/Local/Temp/scp58918/wrk/local2/B22007_nonlinear_scShaper/B20006_Trajectory/antonio_singler.RData")
tbl[,1]

AnnotationScatterPlot(sce,annotation = tbl[,1],dim.reduction.type = "umap",show.legend = TRUE,point.size = 0.3)

