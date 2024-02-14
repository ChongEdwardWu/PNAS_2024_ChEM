### Section 0: Preparation --------------------------------------------------

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

# Load the processed Seurat object
seu <- readRDS("03_KO_pySCNIEC2seurat.rds")

hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), v = T)
hb_genes <- c(grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), v = T))
bad_features <- unique(c(
    hist_genes, hb_genes,
    "Lyz2", # Cre mice
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
        rownames(seu@assays$SCT@counts),
        v = T
    )
))


### Section 1: Data Visulization --------------------------------------------------

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

# Section 2: Annotate celltypes--------------------------------------------------------------

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
        return("B_cell_CD79b")
    } else if (x %in% c(2)) {
        return("SPM_Slamf9")
    } else if (x %in% c(3)) {
        return("LPM_Saa3")
    } else if (x %in% c(4)) {
        return("Neutrophil_Il1b")
    } else if (x %in% c(5)) {
        return("T_cell_Cd3d")
    } else if (x %in% c(6)) {
        return("NK_cell_Nkg7")
    } else if (x %in% c(7)) {
        return("B_cell_Jchain")
    } else if (x %in% c(8)) {
        return("Mast_cell_Cpa3")
    } else {
        return("undefined")
    }
}

# Apply the function to the seurat_clusters vector -l2
seu$CellType_l2 <- factor(sapply(seu$seurat_clusters, assign_cell_type_l2))
table(seu$CellType_l2)

seu$CellType_l2 <- fct_relevel(
    seu$CellType_l2,
   "LPM_Saa3","SPM_Slamf9", "B_cell_CD79b", "B_cell_Jchain", "Neutrophil_Il1b","T_cell_Cd3d","NK_cell_Nkg7","Mast_cell_Cpa3",
)
table(seu$CellType_l2)


# Create a function to assign the cell type based on the value of seurat_clusters -l1
assign_cell_type_l1 <- function(x) {
    if (x %in% c(1,7)) {
        return("B_cell")
    } else if (x %in% c(2, 3)) {
        return("Mf")
    } else if (x %in% c(4)) {
        return("Neutrophil")
    } else if (x %in% c(5)) {
        return("T_cell")
    } else if (x %in% c(6)) {
        return("NK_cell")
    } else if (x %in% c(8)) {
        return("Mast_cell")
    } else {
        return("undefined")
    }
}

# Apply the function to the seurat_clusters vector -l1
seu$CellType_l1 <- factor(sapply(seu$seurat_clusters, assign_cell_type_l1))
levels(seu$CellType_l1)
seu$CellType_l1 <- fct_relevel(seu$CellType_l1,
   "Mf", "B_cell", "Neutrophil","T_cell","NK_cell","Mast_cell",
)
table(seu$CellType_l1)


# Section 3: Identify clusters markers--------------------------------------------------------------
# Load necessary packages for this section
library(future)
plan("multicore", workers = 16)

# Set the default assay and prepare the object for finding markers
Idents(seu) <- "CellType_l1"
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
           file = "results/04_KO_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_DEG", 
           append = F)
jgc()
write.xlsx(seu_ROC %>% as.data.frame(), 
           file = "results/04_KO_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_ROC", 
           append = T)
jgc()
write.xlsx(seu_regulon %>% as.data.frame(), 
           file = "results/04_KO_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_regulon", 
           append = T)
jgc()
write.xlsx(seu_regulon_bin %>% as.data.frame(), 
           file = "results/04_KO_clusters_markers.xlsx", 
           row.names = F,
           sheetName = "clusters_regulon_bin", 
           append = T)
jgc()


### Section 4: Identify differential genes within Mf clusters: KO vs WT  --------------------------------------------------

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

## subset Mf data and prepare
table(seu$CellType_l1)
table(seu$group)
seu_Mf <- seu[,seu$CellType_l1 == "Mf"]

## examine data and prepare
DefaultAssay(seu_Mf) <- "RNA"
seu_Mf$group <- factor(seu_Mf$group)
Idents(seu_Mf) <- "group"
table(seu_Mf$group)
levels(seu_Mf$group)
seu_Mf <- SCTransform(seu_Mf, vst.flavor = "v2")
seu_Mf <- PrepSCTFindMarkers(seu_Mf)


DefaultAssay(seu_Mf) <- "SCT"
features <- setdiff(rownames(seu_Mf@assays$SCT@data), bad_features)

## DEG analysis
resultfile <- "results/04_KOvsWT_Mf-DEG.xlsx"  #! set the result file name here

ChatDEG_SCT <- FindMarkers(
        seu_Mf,
        assay =  "SCT",
        test.use = "wilcox",
        features = features,
        ident.1 = "KO",
        ident.2 = "WT",
        only.pos = FALSE,
        logfc.threshold = 0,
        densify = TRUE
)
write.xlsx(ChatDEG_SCT,
        file = resultfile,
        row.names = T,
        sheetName = "DEG",
        append = FALSE
)
ChatDEG_AUC <- FindMarkers(
        seu_Mf,
        assay =  "AUC",
        test.use = "wilcox",
        ident.1 = "KO",
        ident.2 = "WT",
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
        seu_Mf,
        assay =  "Bin",
        test.use = "wilcox",
        ident.1 = "KO",
        ident.2 = "WT",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0.000,
        densify = TRUE
)
write.xlsx(ChatDEG_Bin,
        file = resultfile,
        row.names = T,
        sheetName = "Binary regulon AUC scores",
        append = TRUE
)

# ! save annotated data
saveRDS(seu_Mf, file = "04_KO_Mfs.rds")
# seu_Mf <- readRDS("04_KO_Mfs.rds")
# Save the Seurat objects
saveRDS(seu, file = "04_KO_all.rds")

### Section 6: Identify differential genes within all clusters: KO vs WT  --------------------------------------------------
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) # For mouse
library(rJava)
library(xlsx)
library(stringr)
library(fgsea)
library(msigdbr)

# Set random seed and multi-threading
set.seed(123)
nworkers <- 16
plan("multicore", workers = nworkers)

# Function to trigger garbage collection in Java and R
jgc <- function() {
    .jcall("java/lang/System", method = "gc")
    gc()
}

# Define cell types to analyze
cell_types <- setdiff(levels(seu$CellType_l1),"Mf")

# Loop through each cell type for analysis
for (cell_type in cell_types) {
    # Set result file name
    resultfile <- paste0("results/04_KOvsWT_", cell_type, "-DEG.xlsx")

    # Extract data for specific cell type
    seu_cell_type <- seu[, seu$CellType_l1 == cell_type]
    DefaultAssay(seu_cell_type) <- "RNA"
    seu_cell_type$group <- factor(seu_cell_type$group)
    Idents(seu_cell_type) <- "group"
    seu_cell_type <- SCTransform(seu_cell_type, vst.flavor = "v2")
    seu_cell_type <- PrepSCTFindMarkers(seu_cell_type)

    # Differential expression analysis
    seu_cell_type$group <- factor(seu_cell_type$group)
    Idents(seu_cell_type) <- "group"
    DefaultAssay(seu_cell_type) <- "SCT"
    features <- setdiff(rownames(seu_cell_type@assays$SCT@data), bad_features)
    ChatDEG_SCT <- FindMarkers(
        seu_cell_type,
        assay = "SCT",
        test.use = "wilcox",
        ident.1 = "KO",
        ident.2 = "WT",
        features = features,
        logfc.threshold = 0,
        densify = TRUE
    )
    ?FindMarkers
    write.xlsx(ChatDEG_SCT,
        file = resultfile,
        row.names = T,
        sheetName = "DEG",
        append = FALSE
    )
    ChatDEG_AUC <- FindMarkers(
        seu_cell_type,
        assay = "AUC",
        test.use = "wilcox",
        ident.1 = "KO",
        ident.2 = "WT",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0,
        densify = TRUE
    )
    write.xlsx(ChatDEG_AUC,
        file = resultfile,
        row.names = T,
        sheetName = "Regulon AUC scores",
        append = TRUE
    )
    ChatDEG_Bin <- FindMarkers(
        seu_cell_type,
        assay = "Bin",
        test.use = "wilcox",
        ident.1 = "KO",
        ident.2 = "WT",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0,
        densify = TRUE
    )
    write.xlsx(ChatDEG_Bin,
        file = resultfile,
        row.names = T,
        sheetName = "Binary regulon AUC scores",
        append = TRUE
    )
    
    # Save analysis results for the current cell type
    save_file_name <- paste0("04_KO_", cell_type, "_DEGresults_", Sys.Date(), ".rds")
    save(ChatDEG_SCT, ChatDEG_AUC, ChatDEG_Bin,
        file = save_file_name
    )

    # Clear variables before next iteration (optional, for memory management)
    rm(ChatDEG_SCT, ChatDEG_AUC, ChatDEG_Bin)

    # Call Java and R garbage collection
    jgc()
}

# Reminder to the user
message("Analysis completed and saved for each cell type.")

# Save the Seurat objects
saveRDS(seu, file = "04_KO_all.rds")
# seu <- readRDS(file = "04_KO_all.rds")
