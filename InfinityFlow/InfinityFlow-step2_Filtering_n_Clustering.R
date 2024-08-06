### Section 0: Preparation --------------------------------------------------

# Clear the environment
rm(list = ls())
graphics.off()
gc()

# # install packages
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("flowWorkspace")
# remotes::install_github("immunedynamics/spectre")
# remotes::install_github("JinmiaoChenLab/Rphenograph")

# BiocManager::install()


# Section 1: Setting up your input data---------------------------------------------------------------

## Load requred packages
source("/ifs1/User/yangbinVM/01.Project/wuchong/.radian_profile")
library(flowCore)
library(flowWorkspace)
library(CATALYST)
library(SingleCellExperiment)
library(stringr)
library(BiocParallel)
library(ggplot2)
library(cowplot)

# Set your working directory
workdir <- "/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/InfinityFlow"
setwd(workdir)
current_date <- format(Sys.Date(), "%y%m%d")
set.seed(123)
nworks <- 32

# Set up input concatenated FACS files derived from InfinityFlow
input_dir <- file.path(workdir, "InfinityFlow_results/")
list.files(input_dir)

# Read in FCS files as flowSet
fcs <- list.files(path = input_dir, pattern = ".fcs$")

fs_all <- read.flowSet(
  fcs,
  path = input_dir,
  transformation = "linearize-with-PnG-scaling",
  truncate_max_range = TRUE
)


# Set sample names and corresponding conditions
sample_names <- sampleNames(fs_all) %>% stringr::str_replace_all(".fcs$", "")
split_names <- stringr::str_split_fixed(sample_names, "_", 3)
condition <- split_names[, 1]
chat_exp <- split_names[, 3]

# Add sample names and conditions to the flowSet object's description
pData(fs_all) <- data.frame(name = sample_names, condition = condition, chat_exp = chat_exp, row.names = sampleNames(fs_all))

# View the number of cells in each experiment
lapply(fs_all@frames, nrow)

# Randomly extract 10000 cells/signals from each experiment, if cells more than 10000
fs <- lapply(fs_all@frames, function(x) {
  n_cells <- nrow(x)
  sample_size <- min(n_cells, 10000)
  x[sample(n_cells, sample_size), ]
})

# Merge the sampled cells/signals into a new flowSet object
fs <- flowSet(frames = fs, phenoData = fs_all@phenoData)


# View the number of cells in each experiment
lapply(fs@frames, nrow)

# Change marker names to facilitate the following process
colnames(fs)
library(stringr)
colnames(fs) <- colnames(fs) %>%
  str_replace_all(" +.*", "") %>%
  str_replace_all("[.-]", "_") %>%
  str_replace_all("\\(", "_aka_") %>%
  str_replace_all("\\)", "") %>%
  str_replace_all("__", "_") %>%
  str_replace_all("GFP_FICT", "ChAT_GFP")

# List unwanted channels
unwantedChn <- c(
  "UMAP1", "UMAP2", "Time", "Exploratory_Ab_ID", "PE_id", # not FACS outcomes
  grep("IgG|_IgM_k", colnames(fs), value = T), # isotypes
  grep("Blank|FJComp", colnames(fs), value = T), # blank channels
  "Zombie_APC_Cy7", "CD45_AF700", "CD45_BV570", "CD45_PC5_5" # barcode
)

# Check which unwanted channels are present in the flowSet
unwantedChn_present <- unwantedChn[unwantedChn %in% colnames(fs)]

# Remove unwanted channels from the flowSet
fs <- fs[, !(colnames(fs) %in% unwantedChn_present)]

# Check the remaining channels
colnames(fs)

# change feature name
grep("VISTA", colnames(fs))
colnames(fs)[grep("VISTA", colnames(fs))[1]] <- "PD_1H_aka_VISTA_1"
colnames(fs)[grep("VISTA", colnames(fs))[2]] <- "PD_1H_aka_VISTA_2"

# Convert the flowSet object to a SingleCellExperiment object
sce <- CATALYST::prepData(fs, FACS = TRUE, transform = TRUE)
# assay(sce, "exprs")[1:5,1:5]

# Create a vector that assigns each cell's sample_name label from each sample
sample_lable <- factor(rep(fs@phenoData$name, times = sapply(fs@frames, nrow)))
# Assign sample_name labels to each cell in the sce object
colData(sce)$sample_id <- factor(sample_lable)

# Create a vector that assigns each cell's group label from each sample
conditions_lable <- factor(rep(condition, times = sapply(fs@frames, nrow)))
# Assign group labels to each cell in the sce object
colData(sce)$condition <- factor(conditions_lable,
  levels = c(
    "Con", "LPS", "Pam"
  )
)


# Create a vector that assigns each cell's chat_exp label from each sample
chat_lable <- factor(rep(chat_exp, times = sapply(fs@frames, nrow)))
# Assign chat_exp labels to each cell in the sce object
colData(sce)$chat_exp <- factor(chat_lable,
  levels = c(
    "Neg", "Pos"
  )
)

# See colData
colData(sce)

# See rownames and drop unwanted parameters
rownames(sce)
head(rowData(sce))
rownames(sce) <- rowData(sce)$marker_name <- rowData(sce)$channel_name

# Add marker information to rowData
rowData(sce)$marker_source <- ifelse(grepl("XGBoost", rownames(sce)), "imputated", "backbone")
table(rowData(sce)$marker_source)
head(rowData(sce))

# Define "type" markers
rowData(sce)$marker_class <- "state"
type_maker_names <- grep("CD11c|CD19_|CD22_|ly_6D|CD20_|MD_1|IgD|IgM|CD11b_PE_CF594|F4_80_APC|MHC_II|Tim_4|CD301a_|Ly_6G|CD49b|CD335|CD127|CD3_|CD3e|TCR_b|TCR_g_d|CD186|FceRIa|CD117", rownames(sce), value = T)
grep("CD11", rownames(sce), value = T)

rowData(sce)$marker_class[rownames(sce) %in% type_maker_names] <- "type"
table(rowData(sce)$marker_class)
# rowData(sce)$marker_class <- "type"

# Rename markers by removing the "_XGBoost_bgc" suffix from the marker names
rownames(sce) <- rowData(sce)$marker_name %>%
  stringr::str_replace_all("_XGBoost_bgc$", "")
rownames(sce)

# Create a new assay named "exprs" that is identical to the "counts" assay in the sce object
# assays(sce)$exprs <- assays(sce)$counts


# Section 2: Cell population identification---------------------------------------------------------------
# run t-SNE/UMAP
library(BiocParallel)
set.seed(123)

# sce <- runDR(sce, "PCA", features = NULL, BPPARAM = MulticoreParam(nworks))
sce <- runDR(sce, "TSNE", features = NULL, BPPARAM = MulticoreParam(nworks)) # features = NULL to use all features
sce <- runDR(sce, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworks)) # features = NULL to use all features



# Section 3: Cell population filtering and annotation - 1st round---------------------------------------------------------------

# Cell population identification with FlowSOM and ConsensusClusterPlus
sce <- cluster(sce,
  features = NULL, # features = NULL to use all features
  xdim = 10, ydim = 10, maxK = 30, seed = 123
)

## ## ## ## ##
# pause here
## ! save sce object
saveRDS(sce, file = "InfinityFlow_Mf-step2_prefiltered_20231010.rds")
# sce <- readRDS("InfinityFlow_Mf-step2_prefiltered_20231010.rds")
## ## ## ## ##

plotExprHeatmap(sce,
  features = "type",
  by = "cluster_id", k = "meta30",
  bars = TRUE, perc = TRUE
)
# manually save the plot as "figures/all_features_Heatmap.png"



# plot the t-SNE and UMAP projections of cells colored by the 20 metaclusters
library(ggplot2)
library(cowplot)

p1 <- plotDR(sce, "TSNE", color_by = "meta30") +
  theme(legend.position = "right")
p2 <- plotDR(sce, "UMAP", color_by = "meta30")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
combined_plot <- plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))

# Define output filenames
png_filename <- "figures/all_features_DR_plot.png"
pdf_filename <- "figures/all_features_DR_plot.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = combined_plot, width = 15, height = 10)

## Cell types annotation
# Create a vector with values 1 to 30
meta30 <- 1:30
# Create a function to assign the cell type based on the value of meta30
assign_cell_type <- function(x) {
  if (x %in% c(1,8,13,15,23)) {
    return("doublets/dead_cells")
  } else if (x %in% c(11,12,17,18,19,20,21,25)) {
    return("B_cells")
  } else if (x %in% c(2,3,4,5,6,7,9,14,16,22)) {
    return("Mfs")
  } else if (x %in% c(10,24)) {
    return("Neutrophils")
  } else if (x %in% c(27,28)) {
    return("abT_cells")
  } else if (x == 30) {
    return("gdT_cells")
  } else if (x %in% c()) {
    return("Unknown")
  } else if (x %in% c(26,29)) {
    return("NK_cells")
  } else {
    return(x)
  }
}

# Apply the function to the meta30 vector
CellType_l1 <- sapply(meta30, assign_cell_type)

# Create the merging_table
merging_table <- data.frame(meta30 = meta30, CellType_l1 = CellType_l1)
merging_table$CellType_l1 <- factor(
  merging_table$CellType_l1
)

# apply manual merging
sce <- mergeClusters(sce,
  k = "meta30",
  table = merging_table,
  id = "CellType_l1",
  overwrite = TRUE
)

# Annotate cells
colData(sce)$CellType_l1 <- cluster_codes(sce)[cluster_ids(sce), "CellType_l1"]
head(colData(sce))

# Comparison of automated and manual merging
plot_grid(
  labels = c("A", "B"),
  plotDR(sce, "UMAP", color_by = "CellType_l1"),
  plotDR(sce, "UMAP", color_by = "meta30")
)

# filter out dead cells or doublets
sce <- sce[, colData(sce)$CellType_l1 != "doublets/dead_cells"]

# Section 4: Cell population filtering and annotation - 2nd round---------------------------------------------------------------
# redo dimension reduction, clustering and annotation
library(BiocParallel)
set.seed(123)
sce <- runDR(sce, "TSNE", features = NULL, BPPARAM = MulticoreParam(nworks)) # features = NULL to use all features
sce <- runDR(sce, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworks)) # features = NULL to use all features

# Cell population identification with FlowSOM and ConsensusClusterPlus
sce <- cluster(sce,
  features = NULL, # features = NULL to use all features
  xdim = 10, ydim = 10, maxK = 20, seed = 123
)

## ## ## ## ##
# pause here
## ! save sce object
saveRDS(sce, file = "InfinityFlow_Mf-step2_1st_filtered_20231010.rds")
# sce <- readRDS("InfinityFlow_Mf-step2_1st_filtered_20231010.rds")
## ## ## ## ##

plotExprHeatmap(sce,
  features = "type",
  by = "cluster_id", k = "meta20",
  bars = TRUE, perc = TRUE
)
# manually save the plot as "figures/filtered_all_features_Heatmap.png"

# plot the t-SNE and UMAP projections of cells colored by the 20 metaclusters
library(ggplot2)
library(cowplot)

# Comparison of automated and manual merging
combined_plot <- plot_grid(
  labels = c("A", "B"),
  plotDR(sce, "UMAP", color_by = "CellType_l1"),
  plotDR(sce, "UMAP", color_by = "meta20")
)

# Define output filenames
png_filename <- "figures/filtered_all_features_DR_plot.png"
pdf_filename <- "figures/filtered_all_features_DR_plot.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = combined_plot, width = 15, height = 10)


# Create a vector with values 1 to 20
meta20 <- 1:20
# Create a function to assign the cell type based on the value of meta20
assign_cell_type_l2 <- function(x) {
if (x %in% c()) {
    return("doublets/dead_cells")
  } else if (x %in% c(3,7,8,13,14)) {
    return("B_cells")
  } else if (x %in% c(9,11,12,15,17,18,19,20)) {
    return("Mfs")
  } else if (x %in% c(5,16)) {
    return("Neutrophils")
  } else if (x %in% c(1)) {
    return("abT_cells")
  } else if (x %in% c(2)) {
    return("gdT_cells")
  } else if (x %in% c(10)) {
    return("DCs")
  } else if (x %in% c(4,6)) {
    return("NK_cells")
  } else {
    return(x)
  }
}

# Apply the function to the meta20 vector
CellType_l2 <- sapply(meta20, assign_cell_type_l2)

# Create the merging_table
merging_table <- data.frame(meta20 = meta20, CellType_l2 = CellType_l2)
merging_table$CellType_l2 <- factor(merging_table$CellType_l2)

# apply manual merging
sce <- mergeClusters(sce,
  k = "meta20",
  table = merging_table,
  id = "CellType_l2",
  overwrite = TRUE
)

# Annotate cells_l2
colData(sce)$CellType_l2 <- cluster_codes(sce)[cluster_ids(sce), "CellType_l2"]
head(cluster_codes(sce))
head(colData(sce))

# filter out dead cells or doublets
sce <- sce[, colData(sce)$CellType_l2 != "doublets/dead_cells"]

# Comparison of automated and manual merging
plot_grid(
  labels = c("A", "B"),
  plotDR(sce, "UMAP", color_by = "meta20"),
  plotDR(sce, "UMAP", color_by = "CellType_l2")
)

# Annotate cells_l1
# Create a function to assign the cell type based on the value of meta15
assign_cell_type_l1 <- function(x) {
  if (grepl("B_cells", x)) {
    return("B_cells")
  } else if (grepl("Mfs", x)) {
    return("Mfs")
  } else if (grepl("NK_cells", x)) {
    return("NK_cells")
  } else if (grepl("Neutrophils", x)) {
    return("Neutrophils")
  } else {
    return(x)
  }
}

CellType_l2 <- levels(colData(sce)$CellType_l2)
# Apply the function to the meta15 vector
CellType_l1 <- sapply(CellType_l2, assign_cell_type_l1)

# Create the merging_table
merging_table <- data.frame(CellType_l2 = CellType_l2, CellType_l1 = CellType_l1)
merging_table$CellType_l1 <- factor(merging_table$CellType_l1)

# apply manual merging
sce <- mergeClusters(sce,
  k = "CellType_l2",
  table = merging_table,
  id = "CellType_l1",
  overwrite = TRUE
)
head(cluster_codes(sce))

# Annotate cells_l1
colData(sce)$CellType_l1 <- cluster_codes(sce)[cluster_ids(sce), "CellType_l1"]


table(colData(sce)$CellType_l1, colData(sce)$CellType_l2)



# redo dimension reduction, clustering and annotation
library(BiocParallel)
library(scran)
library(scater)
set.seed(123)

sce <- fixedPCA(sce,assay.type = "exprs", BPPARAM=MulticoreParam(nworks/4))
reducedDimNames(sce)
sce <- runUMAP(sce, exprs_values ="exprs", dimred="PCA", BPPARAM=MulticoreParam(nworks/4))



sce <- runDR(sce, "TSNE", features = NULL, BPPARAM = MulticoreParam(nworks)) # features = NULL to use all features
sce <- runDR(sce, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworks)) # features = NULL to use all features


# Comparison of automated and manual merging
combined_plot <- plot_grid(
  labels = c("A", "B"),
  plotDR(sce, "UMAP", color_by = "CellType_l1"),
  plotDR(sce, "UMAP", color_by = "CellType_l2")
)

# Define output filenames
png_filename <- "figures/filtered_all_features_DR_plot.png"
pdf_filename <- "figures/filtered_all_features_DR_plot.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = combined_plot, width = 15, height = 10)

# facet by sample id
combined_plot <- plotDR(sce, "UMAP", color_by = "CellType_l1", facet_by = "sample_id")

# Define output filenames
png_filename <- "figures/filtered_all_features_DR_plot_split.png"
pdf_filename <- "figures/filtered_all_features_DR_plot_split.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = combined_plot, width = 15, height = 10)


# ! save data
saveRDS(sce, file = "InfinityFlow_Mf-step2_filtered_20231010.rds")
# sce <- readRDS("InfinityFlow_Mf-step2_filtered_20231010.rds")
