### Section 0: Preparation --------------------------------------------------

suppressMessages(suppressWarnings(source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")))

# # install packages
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }

# BiocManager::install(c(
#   "DropletUtils", "scater", "tidyverse", "scDblFinder", "BiocParallel",
#   "Nebulosa", "ggplot2", "remotes", "robustbase", "future"
# ))

# #remotes::install_github('cole-trapnell-lab/monocle3')
# remotes::install_github(c("satijalab/seurat","satijalab/sctransform","satijalab/seurat-wrappers"), ref = "develop")
# remotes::install_github(c("mojaveazure/seurat-disk"))


# clear the environment
rm(list = ls())
graphics.off()
gc()
source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")

# set your work diretory
workdir <- "/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/R"
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

library(future)
# plan("sequential")
plan("multicore", workers = nworkers)


## ! define samples! ##
samples <- c("WT", "KO", "Pos", "Neg")

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
seurat_ref <- c("WT", "KO")
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
stress_genes <- c("G0s2", "Jun", "Junb", "Jund", "Fos", "Dusp1", "Cdkn1a", "Fosb", "Btg2", "Klf6", "Klf4")
cc_genes <- unique(c(paste0(toupper(substr(tolower(as.character(unlist(cc.genes.updated.2019))), 1, 1)), substr(tolower(as.character(unlist(cc.genes.updated.2019))), 2, nchar(tolower(as.character(unlist(cc.genes.updated.2019)))))))) # , "1700020L24Rik" ,"5730416F02Rik" ,"Agpat4","Asf1b", "Aspm","Ccdc18","Ccr9","Clspn","Cyp4f17","Dek","Dnmt3b" ,"Dtl","Fancm","Fignl1","Gm14214","Gm14730","Gm15428","Gm15448","Gm21992","Gm23248","Gm26465","Gm5145" ,"Mcm4","Mcm8","Mki67","Oip5","Pcna","Pcna-ps2","Pigv","Rnd1","Snrpa","Ube2c"
hist_genes <- grep("Hist", rownames(seu_int@assays$RNA@counts), v = T)
hb_genes <- c(grep("^Hb[ab]-|^HB[^(P)]", rownames(seu_int@assays$RNA@counts), v = T))
bad_features <- unique(c(
    hist_genes, cc_genes, stress_genes, hb_genes,
    "Malat1", "Neat1", "Xist", "Actb",
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
# can skip if already done
ElbowPlot(seu_int, ndims = 50)
# seu_int <- JackStraw(seu_int, num.replicate = 100)
# seu_int <- ScoreJackStraw(seu_int, dims = 1:20)
# JackStrawPlot(seu_int, dims = 1:20)


## UMAP and TSNE Clustering
DefaultAssay(seu_int) <- "integrated"
seu_int <- RunTSNE(seu_int, reduction = "pca", dims = 1:50)
seu_int <- RunUMAP(seu_int, reduction = "pca", dims = 1:50)
seu_int <- FindNeighbors(seu_int, reduction = "pca", dims = 1:50)

## estimate clustering
# try different resolutions
# seu_int@meta.data[,grep("integrated_snn_res.",colnames(seu_int@meta.data),value = T)] <- NULL
seu_int <- FindClusters(
    seu_int,
    resolution = c(seq(0.01, 0.1, 0.01), seq(0.2, 2, 0.1))
    # resolution = c(seq(0.1, 2, 0.1))
)
# relevel cluster names (we do not what cluster "0")
for (i in 1:length(grep("integrated_snn_res.", colnames(seu_int@meta.data)))) {
    j <- grep("integrated_snn_res.", colnames(seu_int@meta.data))[i]
    k <- seu_int@meta.data[, j]
    levels(k) <- as.character(1:length(levels(k)))
    seu_int@meta.data[, j] <- k
}

# estimate clustering in different resolustion
clustree(seu_int@meta.data, prefix = "integrated_snn_res.", return = "plot")

ggsave(
    filename = "figures/01_Clustering_Resolution_all.png",
    width = 250, height = 400, units = "mm",
    dpi = 150, device = "png", bg = "white"
)

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

# ! set group levels
seu_int$group <- factor(seu_int$group,
    levels = c("WT", "KO", "Neg", "Pos")
)

DimPlot(seu_int,
    reduction = "umap", split.by = "group",
    group.by = "seurat_clusters", label = TRUE
) +
    coord_fixed(ratio = 1)

# colnames(seu_int@meta.data)

VlnPlot(seu_int, features = c("sum","detected", "subsets_Mito_percent", "subsets_Rp_percent","subsets_Heatshock_percent"), pt.size = 0)


# find cluster markers
DefaultAssay(seu_int) <- "SCT"
seu_int <- PrepSCTFindMarkers(object = seu_int)
features <- setdiff(rownames(seu_int@assays$SCT@data), bad_features)
# wilcox
seu_DEG <- FindAllMarkers(seu_int,
  only.pos = TRUE,
  min.pct = 0.1,
  features = features,
  logfc.threshold = 0.25,
  densify = TRUE
)
# arrange the results
library(dplyr)
seu_DEG <- seu_DEG %>%
  group_by(cluster) %>%
  arrange(p_val_adj, by_group = T)
# View(seu_DEG)
write.csv(seu_DEG, file = "results/01_seurat_clusters_markers.csv")
DefaultAssay(seu_int) <- "integrated"

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

# ! filter bad cells/clusters, if any, and redo clustering
discaredCl <- c()
seu_int <- seu_int[,!(seu_int$seurat_clusters %in% discaredCl)]


# save the processed seurat object
saveRDS(seu_int,
    file = "01_all_Data_Integration_20230421.rds"
)
# seu_int <- readRDS("01_all_Data_Integration_20230421.rds")
