### Section 0: Preparation --------------------------------------------------

# clear the environment
rm(list = ls())
graphics.off()
gc()

# !set your work diretory
workdir <- "/your_work_dir/data/Peritoneal_ChAT/scRNA_seq"
setwd(workdir)

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
seu_int <- readRDS("01_KO_Data_Integration.rds")

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
loom <- build_loom("KOnWT_int.loom", dgem = exprMat_filter)
# loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)
