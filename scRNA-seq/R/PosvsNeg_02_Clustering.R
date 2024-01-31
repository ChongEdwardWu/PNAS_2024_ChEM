### Section 0: Preparation --------------------------------------------------

# Clearing the R environment and resetting the graphics device
rm(list = ls())
graphics.off()
gc()

# Setting the working directory
workdir <- "/your_work_dir/R"
setwd(workdir)

# Load necessary packages
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(future)
library(rJava)
library(xlsx)
library(stringr)

# Set seed for reproducibility and specify the number of workers for parallel computing
set.seed(123)
nworkers <- 32

# Load the integrated Seurat object
seu <- readRDS("01_WT_Data_Integration.rds")

### Section 1: Re-integrate data --------------------------------------------------
# do this if you filtered cells in the last step: 01_Data_Integration

# Split the Seurat object by 'group'
Seu_list <- SplitObject(seu, split.by = "group")

# Apply SCTransform to each subset using the RNA assay explicitly
Seu_list <- lapply(X = Seu_list, FUN = function(x) {
  x <- SCTransform(object = x, assay = "RNA", verbose = TRUE)
  return(x)
})

# Select features for integration
features <- SelectIntegrationFeatures(object.list = Seu_list, nfeatures = 3000)

# Pre-process the data for integration
Seu_list <- PrepSCTIntegration(
  object.list = Seu_list, 
  anchor.features = features,
  verbose = T)

# Next, identify anchors and integrate the datasets
# !First, set reference, if any
seurat_ref <- c()
# which(names(Seu_list) == seurat.ref)

if (is.null(seurat_ref)) {
    anchors <- FindIntegrationAnchors(
        object.list = Seu_list, normalization.method = "SCT",
        reduction = "rpca", anchor.features = features, verbose = T,
        reference = NULL
    )
} else {
    anchors <- FindIntegrationAnchors(
        object.list = Seu_list, normalization.method = "SCT",
        reduction = "rpca", anchor.features = features, verbose = T,
        reference = which(names(Seu_list) == seurat_ref)
    )
}

seu <- IntegrateData(
    anchorset = anchors, normalization.method = "SCT",
    verbose = T
)

### Section 2: Dimension reduction and clustering --------------------------------------------------
## Run the standard workflow for visualization and clustering
# PCA
hist_genes <- grep("Hist", rownames(seu@assays$RNA@counts), v = T)
hb_genes <- c(grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$RNA@counts), v = T))
bad_features <- unique(c(
    hist_genes, hb_genes,
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
        rownames(seu@assays$RNA@counts),
        v = T
    )
))

PCA_features <-  setdiff(seu@assays$integrated@var.features, bad_features)

seu <- RunPCA(seu,
    npcs = 50, verbose = T,
    features = PCA_features
)

## Harmony integration and UMAP 
set.seed(123)
library(harmony)
library(future)
plan("multicore", workers = nworkers)

seu <- RunHarmony(seu, group.by.vars = "group")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)

## estimate clustering
# try different resolutions
# seu@meta.data[,grep("integrated_snn_res.",colnames(seu@meta.data),value = T)] <- NULL
seu <- FindClusters(
    seu,
    resolution = c(seq(0.01, 0.1, 0.01), seq(0.2, 1, 0.1))
)
# relevel cluster names (we do not what cluster "0")
for (i in 1:length(grep("integrated_snn_res.", colnames(seu@meta.data)))) {
    j <- grep("integrated_snn_res.", colnames(seu@meta.data))[i]
    k <- seu@meta.data[, j]
    levels(k) <- as.character(1:length(levels(k)))
    seu@meta.data[, j] <- k
}

# estimate clustering in different resolustion
library(clustree)
clustree(seu@meta.data, prefix = "integrated_snn_res.", return = "plot")

# ! resolution choice
res <- 0.06

# ! set resolution
Idents(seu) <- seu$seurat_clusters <- factor(seu@meta.data[[paste0("integrated_snn_res.", res)]])
# seu@meta.data[,grep("integrated_snn_res.",colnames(seu@meta.data),value = T)] <- NULL

# Generate a UMAP plot grouped by 'seurat_clusters' and save it to a file
p01 <- DimPlot(seu,
    reduction = "umap", 
    # split.by = "group",
    group.by = "seurat_clusters", label = TRUE
) +
    coord_fixed(ratio = 1)
p01
ggsave(
    p01,
    filename = "figures/02_WT_UMAP.pdf"
)


p02 <- DimPlot(seu,
    reduction = "umap", split.by = "group",
    group.by = "seurat_clusters", label = TRUE
) +
    coord_fixed(ratio = 1)
p02
ggsave(
    p02,
    filename = "figures/02_WT_UMAP_split.pdf"
)

# save the processed seurat object
saveRDS(seu,
    file = "02_WT_Mf_Integration.rds"
)
# seu <- readRDS("02_WT_Mf_Integration.rds")
