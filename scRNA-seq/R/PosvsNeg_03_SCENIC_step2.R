### Section 0: Preparation --------------------------------------------------

# clear the environment
rm(list = ls())
graphics.off()
gc()
source("/your_work_dir/.radian_profile")

# !set your work diretory
workdir <- "/your_work_dir/data/Peritoneal_ChAT/scRNA_seq"
setwd(workdir)

# Section 1: Exploring results in R --------------------------------------------------------------

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)


# !reset your work diretory
SCENICdir <- "/your_work_dir/data/Peritoneal_ChAT/SCENIC/results"
scenicLoomPath <- file.path(SCENICdir, 'seu_int.scenic.loom')
motifEnrichmentFile <- file.path(SCENICdir, 'seu_int.motifs.csv')
regulonAucFile <- file.path(SCENICdir, 'seu_int.auc.csv')
BinarymatFile <- file.path(SCENICdir, 'seu_int.bin.csv')
regulonAucThresholdsFile <- file.path(SCENICdir, 'seu_int.thresholds.csv')

all(file.exists(scenicLoomPath), file.exists(motifEnrichmentFile), file.exists(regulonAucFile),file.exists(BinarymatFile), file.exists(regulonAucThresholdsFile))

# Retrieve AUC scores per cell
loom <- open_loom(scenicLoomPath)
  # Read information from loom file:
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
close_loom(loom)


motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-1,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")
motifEnrichment[1:5,1:5]
length(unique(motifEnrichment$TF))
# View(motifEnrichment)

regulonAUC2  <- data.table::fread(regulonAucFile, header=T, skip=0) %>%
column_to_rownames("Cell") %>%
t()
rownames(regulonAUC2) <- gsub("[(+)]", "", rownames(regulonAUC2))
regulonAUC2[1:5,1:5]
ncol(regulonAUC2)

regulonBin  <- data.table::fread(BinarymatFile, header=T, skip=0) %>%
column_to_rownames("Cell") %>%
t()
rownames(regulonBin) <- gsub("[(+)]", "", rownames(regulonBin))
regulonBin[1:5,1:5]

# get AUC matrix
AUCmat <- AUCell::getAUC(regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))
AUCmat[1:5,1:5]

save(regulons_incidMat, regulons, regulonAUC, motifEnrichment, regulonAUC2, regulonBin,AUCmat,
  file = "03_pySCENIC_results.RData"
)
# load(file = "03_pySCENIC_results.RData")
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Section 2: Incorporate data into seurat object --------------------------------------------------------------

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

packageVersion("SCENIC")

# load seurat object
library('Seurat')
load("03_pySCENIC_results.RData")
seu_int <- readRDS("01_WT_Data_Integration.rds")


# merge AUC matrix into seurat object
seu_int[['AUC']] <- CreateAssayObject(data = regulonAUC2)
DefaultAssay(seu_int) <- 'AUC'
seu_int <- ScaleData(seu_int, assay = 'AUC', features = rownames(regulonAUC2))

# merge AUC-bin matrix into seurat object
seu_int[['Bin']] <- CreateAssayObject(data = regulonBin)
seu_int <- ScaleData(seu_int, assay = 'Bin', features = rownames(regulonBin))

saveRDS(seu_int, file = "03_WT_pySCNIEC2seurat.rds")
