### Section 0: Preparation --------------------------------------------------

# clear the environment
rm(list = ls())
graphics.off()
gc()

suppressMessages(suppressWarnings(source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")))

# # install packages
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }

# BiocManager::install(c(
#   "DropletUtils", "scater", "tidyverse", "scDblFinder", "BiocParallel",
#   "Nebulosa", "ggplot2", "remotes", "robustbase", "future","clustifyr", "clustifyrdatahub"
# ))

# remotes::install_github(c("satijalab/seurat", "satijalab/sctransform","satijalab/seurat-wrappers"), ref = "develop")
# remotes::install_github(c("mojaveazure/seurat-disk", "satijalab/azimuth","ggjlab/scMCA"), ref = "master")

# usethis::create_github_token()
# the current token is: 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
# usethis::edit_r_environ()


### Section 1-identify unqualified clusters ------------------------------------------
suppressMessages(suppressWarnings(library(scater)))

set.seed(123)
nworkers <- 16

## ! type sample name
sample <- "WT"
## ! type sample name
species <- "mm"
# species <- "hs"

## ! path to cellranger outs
cr_path <- file.path("/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/CRcount/WT/outs/filtered_feature_bc_matrix")

## ! path to loom
loom_path <- file.path("/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/CRcount/WT/velocyto/WT.loom")


# set your work diretory
workdir <- file.path("/ifs1/User/yangbinVM/01.Project/wuchong/data/scRNAseq-R", sample)
setwd(workdir)

# load preprocessed data
norm <- readRDS(file = paste0("01_QC_step1_", sample, ".rds"))

# !Summary of filter diagnosis - Discarded clusters
colData(norm)$discard.cl <- ifelse(norm$label %in% c(
  8,21,25,27,29,30,32,34,35,# , # doublets
  14,17,22,36,37 # bad cells
   # RBCs and platelets
),
1, 0
)

# !adjust threholds
qc.nexprs <- norm$detected < 500
qc.mito <- norm$subsets_Mito_percent > 10
qc.rbc <- norm$subsets_RBC_percent > 1
qc.doublet <- norm$DoubletClass == "doublet"
qc.cluster <- norm$discard.cl == 1
discard.new <- qc.nexprs | qc.mito | qc.rbc | qc.doublet | qc.cluster


# Redirect output and messages to the log file
current_date <- format(Sys.Date(), "%Y%m%d")
log_file <- paste0(workdir,"/",sample,"_scRNA_QC_step2_", current_date, ".log")
output_connection <- file(log_file, open = "wt")
sink(output_connection)
sink(output_connection, type = "message")
Sys.time()

# summary of filtering
print("Summary of filtering")
table(discard.new)

norm$discard <- discard.new

# see final filtering results in clusters
sceumap <- gridExtra::grid.arrange(
  plotReducedDim(norm[, norm$discard == 0], "UMAP", colour_by = "label", text_by = "label") +
    ggtitle("Cells remained"),
  plotReducedDim(norm[, norm$discard == 1], "UMAP", colour_by = "label", text_by = "label") +
    ggtitle("Cells discarded"),
  ncol = 2
)
ggsave(
  filename = "figures/08_Final_Filtering_UMAP.png",
  plot = sceumap,
  width = 250, height = 250, units = "mm",
  dpi = 150, device = "png", bg = "white"
)

# We can diagnose cell type loss by looking for systematic
# differences in gene expression between the discarded and retained cells.
rownames(rowData(norm)) <- rowData(norm)$Symbol
lost <- calculateAverage(counts(norm)[, norm$discard])
kept <- calculateAverage(counts(norm)[, !norm$discard])
suppressMessages(suppressWarnings(library(edgeR)))
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)
png(file = "figures/09_filtered_cells_DEGs.png", width = 250, height = 250, units = "mm", res = 150)
plot(abundance, logFC, xlab = "Average count", ylab = "Log-FC (lost/kept)", pch = 16)
points(abundance[discard.new], logFC[discard.new], col = "dodgerblue", pch = 16)
dev.off()

# The presence of a distinct population in the discarded pool
# as a set of genes that are strongly upregulated in lost
# Let's see what these genes are
DEGs <- as.data.frame(lost[logFC > 1])
# View(DEGs)
print("DEGs between cells remained and discarded")
rowData(norm)$Symbol[which(rowData(norm)$ID %in% rownames(DEGs))]

### Section 2-final cell filtering decision ------------------------------------------

## ## ## ## ## ## ## ## ## ## ## ## #
##  Final Cell Filtering Decision ##
## ## ## ## ## ## ## ## ## ## ## ## #
suppressMessages(suppressWarnings(library(scran)))
suppressMessages(suppressWarnings(library(BiocParallel)))

# Re-subset our SingleCellExperiment object to retain only the detected cells
set.seed(123)

# Filter cells
sce <- norm[, (!norm$discard)]

# Assign sample
sce$group <- sample

# Assign Cell-Cycle Scores
if (species == "mm") {
  CCpairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package = "scran"))
} else if(species == "hs") {
  CCpairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
}

assignments <- cyclone(sce, CCpairs,
  gene.names = rowData(sce)$ID,
  BPPARAM = MulticoreParam(nworkers)
)
sce$G1.score <- assignments$normalized.scores$G1
sce$S.score <- assignments$normalized.scores$S
sce$G2M.score <- assignments$normalized.scores$G2M
sce$CCphase <- assignments$phases

# calculate difference between the G2M and S phase scores
sce$CC.Diff <- sce$S.score - sce$G2M.score
# Add integrated cell IDs
library(stringr)
sce$CellID <- str_c(paste0(sample, ":"), str_sub(string = colnames(sce), start = 1, end = 16), "x")
colnames(sce) <- sce$CellID

## Add cell type annotation -- For mouse only
if (species == "mm") {
  suppressMessages(suppressWarnings(library(scMCA)))
  rownames(sce) <- rowData(sce)$Symbol
  mca_result <- scMCA(
    scdata = assays(sce)@listData[["logcounts"]],
    numbers_plot = 3
  )
  # scMCA_vis(mca_result)
  # View(mca_result[["scMCA_probility"]])
  sce$CellType_mca <- mca_result[["scMCA"]]

  # perform annotation by calling SingleR()
  suppressMessages(suppressWarnings(library(SingleR)))
  load("/ifs1/User/yangbinVM/01.Project/wuchong/Ref_genome/CellType_ref/Immgen_ref/ImmGen_reference-Heng_2008.RData")
  pred <- SingleR(
    test = sce,
    ref = immgen,
    labels = immgen$label.fine,
    assay.type.test = 1,
    BPPARAM = MulticoreParam(nworkers)
  )
  sce$CellType_immgen <- mca_result[["scMCA"]] <- pred$pruned.labels
}

## export colDatas
colData <- as_tibble(sce@colData)

# save colData
saveRDS(
  colData,
  file = paste0(current_date, "_", sample, "_colData.rds")
)
# colData <- readRDS(paste0(current_date,"_",sample,"_colData.rds"))

# Section 3: construct Seurat object ----------------------------------------------
## Using Seurat to integrate datasets ##
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SeuratDisk)))
suppressMessages(suppressWarnings(library(SeuratWrappers)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(ggplot2)))

library(future)
# plan('sequential')
plan('multicore', workers = nworkers)

set.seed(123)

# read in 10x cellranger outcomes
counts <- Read10X(data.dir = cr_path)
seu <- CreateSeuratObject(counts = counts)

# rename cells in accordance with the loom file
newnames <- str_c(paste0(sample, ":"), str_sub(string = Cells(seu), start = 1, end = 16), "x")
seu <- RenameCells(seu, new.names = newnames)
print("Check whether the Cell IDs are properly formatted in the created seurat object")
all(Cells(seu) %in% colData$CellID)
# filter cells
seu <- seu[, Cells(seu) %in% colData$CellID]

# read in loom file
if (file.exists(loom_path)) {
  loom <- ReadVelocity(file = loom_path)
  # convert loom to seurat object
  seu_loom <- as.Seurat(loom)
  seu_loom <- subset(seu_loom, cells= Cells(seu))
  # create spliced and unsliced assays
  spliced_assay <- CreateAssayObject(counts = seu_loom@assays$spliced$counts[, Cells(seu)])
  unspliced_assay <- CreateAssayObject(counts = seu_loom@assays$unspliced$counts[, Cells(seu)])
  print("Check whether the Cell IDs are the same in the seurat object and in the newly created assays")
  all(Cells(seu) == colnames(spliced_assay), Cells(seu) == colnames(unspliced_assay))
  seu[["spliced"]] <- spliced_assay
  seu[["unspliced"]] <- unspliced_assay
}

# import the colData from the above SingleCellExperiment subject
seu$CellID <- Cells(seu)
seu@meta.data <- left_join(
  as_tibble(seu@meta.data),
  colData,
  by = "CellID"
) %>%
  column_to_rownames("CellID")

## Perform normalization and dimensionality reduction
# perform normalization using SCTransform
# with an additional flag vst.flavor="v2" to invoke the v2 regularization
vars_to_regress <- c("subsets_Mito_percent", "subsets_Rp_percent", "S.score", "G2M.score") # total.UMI


seu <- SCTransform(seu,
  vst.flavor = "v2",
  verbose = T,
  return.only.var.genes = F,
  # We suggest regressing out the difference between the G2M and S phase scores.
  # This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase among proliferating cells (which are often uninteresting), will be regressed out of the data
  vars.to.regress = vars_to_regress
) %>% RunPCA(npcs = 30, verbose = T)


# Cell type annotation - Azimuth
library(Azimuth)
options(timeout=10000)
if (species == "hs") {
  # Map data to Azimuth annotation
  # The RunAzimuth function can take a Seurat object as input
  seu.annot <- RunAzimuth(
    seu,
    reference = "/ifs1/User/yangbinVM/01.Project/wuchong/Ref_genome/CellType_ref/Azimuth_ref/Human_PBMC"
  )
  # transfer the annotation results (which is in the metadata)
  seu@meta.data <- seu.annot@meta.data
}

# save result and print sucess information
saveRDS(
  seu,
  file = paste0(sample, "_seu.rds")
)
print(paste0("Job ",sample," has been sucessfully done at ", Sys.time(), "."))

sink()
sink(type = "message")

rm(list = ls())
graphics.off()
gc()

# Run in server
# nohup Rscript --vanilla /ifs1/User/yangbinVM/01.Project/wuchong/data/scRNAseq-R/WT/WT_scRNA_QC_step2.R &

