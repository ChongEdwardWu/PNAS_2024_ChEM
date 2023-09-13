### Section 0: Preparation --------------------------------------------------

suppressMessages(suppressWarnings(source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")))

# # install packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::version()

# ## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3","network")) # Optional. Can be replaced by GRNBoost

# ## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# BiocManager::install(c("viridis","ComplexHeatmap"))
# # To export/visualize in http://scope.aertslab.org
# remotes::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# # You are now ready to install SCENIC:
# remotes::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")

# clear the environment
rm(list = ls())
graphics.off()
gc()
source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")

# !set your work diretory
workdir <- "/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/R"
setwd(workdir)



# if Github token is inactivated
# usethis::create_github_token()
# the current token is: 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
# usethis::edit_r_environ()


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

packageVersion("SCENIC")


# !reset your work diretory
SCENICdir <- "/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/SCENIC/results"
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
  file = "03_pySCENIC_results_20230430.RData"
)
# load(file = "03_pySCENIC_results_20230427.RData")
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
load("03_pySCENIC_results_20230430.RData")
seu_int <- readRDS("01_all_Data_Integration_20230421.rds")


# merge AUC matrix into seurat object
seu_int[['AUC']] <- CreateAssayObject(data = regulonAUC2)
DefaultAssay(seu_int) <- 'AUC'
seu_int <- ScaleData(seu_int, assay = 'AUC', features = rownames(regulonAUC2))

# merge AUC-bin matrix into seurat object
seu_int[['Bin']] <- CreateAssayObject(data = regulonBin)
seu_int <- ScaleData(seu_int, assay = 'Bin', features = rownames(regulonBin))

saveRDS(seu_int, file = "01.5_all_Clustering_20230430.rds")
