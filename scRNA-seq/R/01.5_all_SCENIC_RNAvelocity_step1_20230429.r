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


# Section 1: Load Seurat Obj---------------------------------------------------------------

# load libraries
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)
set.seed(123)

# !load integrated seurat object
seu_int <- readRDS("01_all_Data_Integration_20230421.rds")
seu_int <- SCTransform(seu_int, vst.flavor = "v2")

# Section 2: Filtering genes for calculating time---------------------------------------------------------------
exprMat <- seu_int@assays$SCT@data
cellInfo <- data.frame(seu_int@meta.data)
dim(cellInfo)
colnames(cellInfo)

loci1 <- which(rowSums(exprMat) > log1p(1 * .01 * ncol(exprMat)))
table(rowSums(exprMat) > log1p(1 * .01 * ncol(exprMat)))
exprMat_filter <- exprMat[loci1, ]
dim(exprMat_filter)

# Section 3: Tranfer into loom files---------------------------------------------------------------
loom <- build_loom("seu_int.loom", dgem = exprMat_filter)
# loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)



### Section 4: convert Seurat object to h5ad --------------------------------------------------
## load packages
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
set.seed(123)

# integrated
DefaultAssay(seu_int) <- "RNA"
SaveH5Seurat(seu_int, filename = "Int.h5Seurat", overwrite = T)
Convert("Int.h5Seurat", dest = "h5ad", overwrite = T)

# PreCAR
PreCAR <- seu_int[,seu_int$group == "PreCAR"]
PreCAR
DefaultAssay(PreCAR) <- "RNA"
SaveH5Seurat(Pos, filename = "PreCAR.h5Seurat", overwrite = T)
Convert("PreCAR.h5Seurat", dest = "h5ad", overwrite = T)

# PreCC5
PreCC5 <- seu_int[,seu_int$group == "PreCC5"]
PreCC5
DefaultAssay(PreCC5) <- "RNA"
SaveH5Seurat(PreCC5, filename = "PreCC5.h5Seurat", overwrite = T)
Convert("PreCC5.h5Seurat", dest = "h5ad", overwrite = T)



# PostCAR
PostCAR <- seu_int[,seu_int$group == "PostCAR"]
PostCAR
DefaultAssay(PostCAR) <- "RNA"
SaveH5Seurat(PostCAR, filename = "PostCAR.h5Seurat", overwrite = T)
Convert("PostCAR.h5Seurat", dest = "h5ad", overwrite = T)

# PostCC5
PostCC5 <- seu_int[,seu_int$group == "PostCC5"]
PostCC5
DefaultAssay(PostCC5) <- "RNA"
SaveH5Seurat(PostCC5, filename = "PostCC5.h5Seurat", overwrite = T)
Convert("PostCC5.h5Seurat", dest = "h5ad", overwrite = T)


