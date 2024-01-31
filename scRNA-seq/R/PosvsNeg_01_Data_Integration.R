### Section 0: Preparation --------------------------------------------------

suppressMessages(suppressWarnings(source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")))

# set your work diretory
workdir <- "/your_work_dir/R"
setwd(workdir)

# creat diretories for figures and results
if (file.exists(file.path(workdir, "figures"))) {
} else {
    dir.create(file.path(workdir, "figures"))
}
if (file.exists(file.path(workdir, "results"))) {
} else {
    dir.create(file.path(workdir, "results"))
}

### Section 1: Data integration --------------------------------------------------

## Using Seurat to integrate datasets ##
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(clustree)
set.seed(123)
nworkers <- 32

## ! define samples! ##
samples <- c("WT", "Pos", "Neg")

## Define seurat list
Seu_list <- list()

for (i in 1:length(samples)) {
    seu <- readRDS(file = file.path(workdir,"samples", samples[i], paste0(samples[i], "_seu.rds")))
    Seu_list <- append(Seu_list, seu)
    rm(seu)
}
names(Seu_list) <- samples

# Next, select features for downstream integration
features <- SelectIntegrationFeatures(object.list = Seu_list, nfeatures = 3000)

# Run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
Seu_list <- PrepSCTIntegration(
    object.list = Seu_list, anchor.features = features,
    verbose = T
)

# Next, identify anchors and integrate the datasets
# !First, set reference, if any
seurat_ref <- c("WT")
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

seu_int <- IntegrateData(
    anchorset = anchors, normalization.method = "SCT",
    verbose = T
)


### Section 2: Dimension reduction and clustering --------------------------------------------------
## Run the standard workflow for visualization and clustering

DefaultAssay(seu_int) <- "integrated"

# PCA
hist_genes <- grep("Hist", rownames(seu_int@assays$RNA@counts), v = T)
hb_genes <- c(grep("^Hb[ab]-|^HB[^(P)]", rownames(seu_int@assays$RNA@counts), v = T))
bad_features <- unique(c(
    hist_genes, hb_genes,
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
        rownames(seu_int@assays$RNA@counts),
        v = T
    )
))

PCA_features <-  setdiff(seu_int@assays$integrated@var.features, bad_features)

seu_int <- RunPCA(seu_int,
    npcs = 50, verbose = T,
    features = PCA_features
)

# Determine the 'dimensionality' of the dataset
ElbowPlot(seu_int, ndims = 50)

## UMAP and TSNE Clustering
library(future)
# plan("sequential")
plan("multicore", workers = nworkers)
set.seed(123)
DefaultAssay(seu_int) <- "integrated"
seu_int <- RunTSNE(seu_int, reduction = "pca", dims = 1:30)
seu_int <- RunUMAP(seu_int, reduction = "pca", dims = 1:30)
seu_int <- FindNeighbors(seu_int, reduction = "pca", dims = 1:30)

## estimate clustering
# try different resolutions
seu_int <- FindClusters(
    seu_int,
    resolution = c(seq(0.01, 0.1, 0.01), seq(0.2, 1, 0.1))
)
# relevel cluster names (we do not what cluster "0")
for (i in 1:length(grep("integrated_snn_res.", colnames(seu_int@meta.data)))) {
    j <- grep("integrated_snn_res.", colnames(seu_int@meta.data))[i]
    k <- seu_int@meta.data[, j]
    levels(k) <- as.character(1:length(levels(k)))
    seu_int@meta.data[, j] <- k
}

# estimate clustering in different resolustion
library(clustree)
clustree(seu_int@meta.data, prefix = "integrated_snn_res.", return = "plot")

# ! resolution choice
res <- 0.3

# Plot 1: see group difference in clusters-UMAP
DimPlot(seu_int,
    reduction = "umap", split.by = "group",
    group.by = paste0("integrated_snn_res.", res), label = TRUE
) +
    coord_fixed(ratio = 1)

# ! set resolution
Idents(seu_int) <- seu_int$seurat_clusters <- factor(seu_int@meta.data[[paste0("integrated_snn_res.", res)]])
# seu_int@meta.data[,grep("integrated_snn_res.",colnames(seu_int@meta.data),value = T)] <- NULL

# check cell QC
VlnPlot(seu_int, features = c("sum","detected", "subsets_Mito_percent", "subsets_Rp_percent","subsets_Heatshock_percent"), pt.size = 0)

# View cluster annotation
cluster_annot <- tibble(
    seu_int$seurat_clusters,
    seu_int$group,
    seu_int$CellType_immgen
) %>%
    set_names("Cluster", "Source", "CellType") %>%
    group_by(Cluster, CellType) %>%
    summarise(no.cell = n()) %>%
    group_by(Cluster) %>%
    mutate(
        total.no = sum(no.cell),
        perc = 100 * no.cell / total.no
    ) %>%
    arrange(Cluster, dplyr::desc(perc)) %>%
    top_n(n = 5, wt = perc)

View(cluster_annot)

# View and plot cluster distribution
source_cluster <- tibble(
    seu_int$seurat_clusters,
    seu_int$group
) %>%
    set_names("Cluster", "Group") %>%
    group_by(Cluster, Group) %>%
    summarise(no.cell = n()) %>%
    group_by(Group) %>%
    mutate(
        total.no = sum(no.cell),
        perc = 100 * no.cell / total.no
    ) %>%
    dplyr::select(Cluster, Group, perc)
View(source_cluster %>%
    spread(Cluster, perc))

library(ggplot2)
ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
    geom_col(colour = "black") +
    # scale_fill_manual(values=cl) +
    coord_fixed(ratio = 1 / 10) +
    theme_bw() +
    xlab("group") +
    ylab("%")


### Section 3: filter bad cells/clusters and or non-macrophages, if any --------------------------------------------------

# ! filter bad cells/clusters and or non-macrophages, if any
discardCl <- seu_int$seurat_clusters %in% c(2,7,8,10,12,13)
table(discardCl)
discardDim <- ((seu_int[["umap"]]@cell.embeddings[,1] > 2) & (seu_int[["umap"]]@cell.embeddings[,2] < -2) |(seu_int[["umap"]]@cell.embeddings[,2] < -6))
table(discardDim)

table((discardCl | discardDim))
seu_filt <- seu_int[,!(discardCl | discardDim)]

DimPlot(seu_filt,
    reduction = "umap", split.by = "group",
    group.by = "seurat_clusters", label = TRUE
) +
    coord_fixed(ratio = 1)


# save the processed seurat object
saveRDS(seu_filt,
    file = "01_WT_Data_Integration.rds"
)
# seu_filt <- readRDS("01_WT_Data_Integration.rds")
