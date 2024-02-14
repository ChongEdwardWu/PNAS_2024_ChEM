### Section 0: Preparation --------------------------------------------------

# Suppressing messages and warnings when sourcing the radian profile script
suppressMessages(suppressWarnings(source("/your_work_dir/.radian_profile")))

# Clearing the R environment and resetting the graphics device
rm(list = ls())
graphics.off()
gc()

# Setting the working directory
workdir <- "/your_work_dir/data/Peritoneal_ChAT/scRNA_seq"
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
library(future)
# plan("sequential")
plan("multicore", workers = nworkers)

# Load a previously saved Seurat object for further analysis
seu <- readRDS("03_WT_pySCNIEC2seurat.rds")

hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), v = T)
hb_genes <- c(grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), v = T))
bad_features <- unique(c(
    hist_genes, hb_genes,
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
        rownames(seu@assays$SCT@counts),
        v = T
    )
))

### Section 1: Data Visualization --------------------------------------------------

# Generate a UMAP plot grouped by 'seurat_clusters' and save it to a file
p01 <- DimPlot(seu,
    reduction = "umap", 
    # split.by = "group",
    group.by = "seurat_clusters", label = TRUE
) +
    coord_fixed(ratio = 1)
p01


p02 <- DimPlot(seu,
    reduction = "umap", split.by = "group",
    group.by = "seurat_clusters", label = TRUE
) +
    coord_fixed(ratio = 1)
p02

# Section 2: Identify clusters markers--------------------------------------------------------------
# Load necessary packages for this section
library(future)
plan("multicore", workers = 32)

# Set the default assay and prepare the object for finding markers
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "SCT"
seu <- PrepSCTFindMarkers(object = seu)

features <- setdiff(rownames(seu@assays$SCT@data), bad_features)

# wilcox
seu_DEG <- FindAllMarkers(seu,
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

# ROC
seu_ROC <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    features = features,
    test.use = "roc",
    densify = TRUE
)

# regulon
seu_regulon <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    features = NULL,
    test.use = "roc",
    assay =  "AUC",
    logfc.threshold = 0.005,
    densify = TRUE
)

# regulon_Bin
seu_regulon_bin <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    features = NULL,
    test.use = "roc",
    assay =  "Bin",
    logfc.threshold = 0.005,
    densify = TRUE
)

# for exporting results
library(rJava)
library(xlsx)
library(stringr)
jgc <- function()
{
  .jcall("java/lang/System", method ="gc")
  gc()
  
}
write.xlsx(seu_DEG %>% as.data.frame(), 
           file = "results/04_WT_Mf_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_DEG", 
           append = F)
jgc()
write.xlsx(seu_ROC %>% as.data.frame(), 
           file = "results/04_WT_Mf_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_ROC", 
           append = T)
jgc()
write.xlsx(seu_regulon %>% as.data.frame(), 
           file = "results/04_WT_Mf_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_regulon", 
           append = T)
jgc()
write.xlsx(seu_regulon_bin %>% as.data.frame(), 
           file = "results/04_WT_Mf_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_regulon_bin", 
           append = T)
jgc()

# Section 4: Annotate Cell Types --------------------------------------------------------------
# Annotate clusters with cell types based on specific markers
# View cluster annotation
cluster_annot <- tibble(
    seu$seurat_clusters,
    seu$group,
    seu$CellType_mca
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


# Create a function to assign the cell type based on the value of seurat_clusters -l2
assign_cell_type_l2 <- function(x) {
  if (x %in% c(1)) {
    return("SPM_Itgal")
  } else if (x %in% c(2)) {
    return("LPM_Gata6")
  } else if (x %in% c(3)) {
    return("SPM_Lyz1")
  } else if (x %in% c(4)) {
    return("SPM_Cd74")
  } else if (x %in% c(5)) {
    return("SPM_Mki67")
  } else {
    return("undefined")
  }
}

# Apply the function to the seurat_clusters vector -l2
seu$CellType_l2 <- factor(sapply(seu$seurat_clusters, assign_cell_type_l2))
table(seu$CellType_l2)

seu$CellType_l2 <- fct_relevel(
    seu$CellType_l2,
    "LPM_Gata6", "SPM_Itgal", "SPM_Lyz1", "SPM_Cd74", "SPM_Mki67"
)
table(seu$CellType_l2)


# Create a function to assign the cell type based on the value of seurat_clusters -l1
assign_cell_type_l1 <- function(x) {
  if (grepl("LPM_",x)) {
    return("LPM")
  } else if (grepl("SPM_",x)) {
    return("SPM")
  }
}

# Apply the function to the seurat_clusters vector -l1
seu$CellType_l1 <- factor(sapply(seu$CellType_l2, assign_cell_type_l1))
levels(seu$CellType_l1)
seu$CellType_l1 <- fct_relevel(seu$CellType_l1,
    "LPM","SPM")
table(seu$CellType_l1)


### Section 5: Identify differential genes within SPM clusters: Pos vs Neg  --------------------------------------------------

## load packages and prepare environment
library(clusterProfiler)
#library(org.Hs.eg.db) # for human
library(org.Mm.eg.db) # for mouse
library(rJava)
library(xlsx)
library(stringr)

jgc <- function(){
  .jcall("java/lang/System", method ="gc")
  gc()
  
}
plan("multicore", workers = nworkers)
set.seed(123)

## subset SPM data and prepare
table(seu$CellType_l1)
seu_SPM <- seu[,seu$CellType_l1 == "SPM"]

## examine data and prepare
DefaultAssay(seu_SPM) <- "RNA"
seu_SPM$group <- factor(seu_SPM$group)
Idents(seu_SPM) <- "group"
table(seu_SPM$group)
levels(seu_SPM$group)
seu_SPM <- SCTransform(seu_SPM, vst.flavor = "v2")
seu_SPM <- PrepSCTFindMarkers(seu_SPM)


DefaultAssay(seu_SPM) <- "SCT"
features <- setdiff(rownames(seu_SPM@assays$SCT@data), bad_features)

## DEG analysis
resultfile <- "results/04_PosvsNeg_SPM-DEG.xlsx"  #! set the result file name here

ChatDEG_SCT <- FindMarkers(
        seu_SPM,
        assay =  "SCT",
        test.use = "wilcox",
        ident.1 = "Pos",
        features = features,
        ident.2 = "Neg",
        min.pct = 0.01,
        only.pos = FALSE,
        logfc.threshold = 0,
        recorrect_umi = FALSE,
        densify = F
)
write.xlsx(ChatDEG_SCT,
        file = resultfile,
        row.names = T,
        sheetName = "DEG",
        append = FALSE
)
ChatDEG_AUC <- FindMarkers(
        seu_SPM,
        assay =  "AUC",
        test.use = "wilcox",
        ident.1 = "Pos",
        ident.2 = "Neg",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0.000,
        densify = TRUE
)
write.xlsx(ChatDEG_AUC,
        file = resultfile,
        row.names = T,
        sheetName = "Regulon AUC scores",
        append = TRUE
)
ChatDEG_Bin <- FindMarkers(
        seu_SPM,
        assay =  "Bin",
        test.use = "wilcox",
        ident.1 = "Pos",
        ident.2 = "Neg",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0.005,
        densify = TRUE
)
write.xlsx(ChatDEG_Bin,
        file = resultfile,
        row.names = T,
        sheetName = "Binary regulon AUC scores",
        append = TRUE
)



# ! save annotated data
saveRDS(seu, file = "04_WT_Mfs.rds")
saveRDS(seu_SPM, file = "04_WT_SPMs.rds")
