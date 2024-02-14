### Section 0: Preparation --------------------------------------------------

# Clear the environment
rm(list = ls())
graphics.off()
gc()

# Section 1: Setting up your input data---------------------------------------------------------------

## Load requred packages
library(flowCore)
library(flowWorkspace)
library(CATALYST)
library(SingleCellExperiment)
library(stringr)
library(BiocParallel)
library(ggplot2)
library(cowplot)

# Set your working directory
workdir <- "/your_work_dir/data/Peritoneal_ChAT/InfinityFlow"
setwd(workdir)
current_date <- format(Sys.Date(), "%y%m%d")
set.seed(123)
nworkers <- 16

# ! load data
sce <- readRDS("InfinityFlow_Mf-step2_filtered.rds")

# construct colData
head(colData(sce))
colData(sce)$CellType <- factor(colData(sce)$CellType_l2)
colData(sce)$ChAT_CellType <- factor(paste0(colData(sce)$chat_exp, "-", colData(sce)$CellType))
colData(sce)$txt_ChAT_CellType <- factor(paste0(colData(sce)$condition, "-", colData(sce)$ChAT_CellType))
head(colData(sce))
levels(colData(sce)$txt_ChAT_CellType)
table(colData(sce)$condition, colData(sce)$ChAT_CellType)

levels(colData(sce)$CellType)

# Section 2: Find cell type markers ---------------------------------------------------------------
library(scran)
library(scater)
# for writing results
library(rJava)
library(xlsx)
library(stringr)
jgc <- function() {
    .jcall("java/lang/System", method = "gc")
    gc()
}

marker_info <- scoreMarkers(sce, colData(sce)$CellType, full.stats = FALSE, assay.type = "exprs")
marker_info
selected_cols <- c("mean.AUC", "min.AUC", "mean.logFC.cohen", "min.logFC.cohen", "self.average", "other.average")

# Get cell types
cell_types <- levels(colData(sce)$CellType)

# Initialize 'append' parameter
append <- FALSE

cell_type_mks <- c()
# Loop through each cell type
for (cell_type in cell_types) {
  
  # Check if the cell_type exists in marker_info
  if (cell_type %in% names(marker_info)) {
    
    # Extract and process the data for the current cell type
    cell_data <- marker_info[[cell_type]] %>% 
      as.data.frame() %>% 
      arrange(desc(mean.AUC)) %>% 
      mutate(discriminating_marker = ifelse(min.AUC > 0.7 & mean.logFC.cohen > 1, "Yes", "No")) %>% 
      dplyr::select(all_of(selected_cols), discriminating_marker)
    
    # Write data to Excel, with 'append' parameter
    write.xlsx(cell_data,
               file = "results/CellType_markers.xlsx",
               row.names = TRUE,
               sheetName = str_replace_all(cell_type, pattern = "[^[:alnum:]]+", replacement = "_"),
               append = append
    )
    
    # Ensure that subsequent sheets are appended
    append <- TRUE
    .jcall("java/lang/System", method = "gc")
    gc()

    # extract top markers
    cell_type_mks <- c(cell_type_mks,
    cell_data %>%
      dplyr::arrange(desc(min.AUC)) %>%
      dplyr::filter(min.AUC > 0.7) %>%
      dplyr::filter(mean.logFC.cohen > 1) %>%
      rownames())

  }
}

cell_type_mks <- unique(cell_type_mks)

# plot the cell type markers heatmap
library(ggplot2)

p00 <- plotExprHeatmap(sce,
  features = cell_type_mks,
  by = "cluster_id", 
  k = "CellType_l1",
  scale = "last",
  bars = TRUE, 
  perc = TRUE,
  hm_pal = rev(hcl.colors(10, "YlGnBu"))
)
p00

# Define output filenames
png_filename <- "figures/FigS1_cell_type_markers_heatmap.png"
pdf_filename <- "figures/FigS1_cell_type_markers_heatmap.pdf"

# !Save the combined plot manually

# Section 3: DR plots ---------------------------------------------------------------

# Define "type" markers
rowData(sce)$marker_class <- "state"
rowData(sce)$marker_class[rownames(sce) %in% cell_type_mks] <- "type"
table(rowData(sce)$marker_class)
# rowData(sce)$marker_class <- "type"

# run t-SNE/UMAP using "type" features only
library(BiocParallel)
set.seed(123)
sce <- runDR(sce, "TSNE", features = "type", BPPARAM = MulticoreParam(nworkers), name = "TSNE_type")
sce <- runDR(sce, "UMAP", features = "type", BPPARAM = MulticoreParam(nworkers), name = "UMAP_type")


# plot the t-SNE or UMAP projections of cells
library(ggplot2)
library(cowplot)

# T-SNE
p01 <- plot_grid(
  labels = c("A", "B"),
  plotDR(sce, "TSNE_type", color_by = "CellType_l1"),
  plotDR(sce, "TSNE_type", color_by = "meta20")
)
p01


# Define output filenames
png_filename <- "figures/Fig1a_TSNE_plot.png"
pdf_filename <- "figures/Fig1a_TSNE_plot.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p01, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p01, width = 15, height = 10)


p02 <- plot_grid(
  labels = c("A", "B"),
  plotDR(sce, "UMAP_type", color_by = "CellType_l1"),
  plotDR(sce, "UMAP_type", color_by = "meta20")
)
p02
# Define output filenames
png_filename <- "figures/Fig1b_UMAP_plot.png"
pdf_filename <- "figures/Fig1b_UMAP_plot.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p02, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p02, width = 15, height = 10)


# Section 4: DM plot - facet by sample id ---------------------------------------------------------------
colData(sce)$sample_id <- factor(colData(sce)$sample_id,
  levels = c("Con_GFP_Pos","LPS_GFP_Pos", "Pam_GFP_Pos", "Con_GFP_Neg", "LPS_GFP_Neg",  "Pam_GFP_Neg")
)
p03 <- plotDR(sce, "TSNE_type", color_by = "CellType", facet_by = "sample_id")
# Define output filenames
png_filename <- "figures/Fig2a_TSNE_plot_split.png"
pdf_filename <- "figures/Fig2a_TSNE_plot_split.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p03, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p03, width = 15, height = 10)

# 2D hexagonal histogram plot
# Extract t-SNE coordinates
tsne_coords <- reducedDim(sce, "TSNE_type")
# Convert to a dataframe
tsne_df <- as.data.frame(tsne_coords)
# Add sample information
tsne_df$sample_id <- sce$sample_id
# Create a 2D hexagonal histogram plot
p04 <- ggplot(tsne_df, aes(x = TSNE1, y = TSNE2)) +
  geom_hex(bins = 110) +                    # Use hexagonal bins with a specified number
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 25)
    ) +                 # Use Viridis color scale
  facet_wrap(~sample_id) +                 # Facet plot by sample_id
  theme_bw() +                             # Use a white background theme
  theme(panel.spacing = unit(1, "lines"))  # Adjust spacing between panels
p04

# Define output filenames
png_filename <- "figures/Fig2b_TSNE_plot_density.png"   # Filename for PNG output
pdf_filename <- "figures/Fig2b_TSNE_plot_density.pdf"   # Filename for PDF output
# Save the plot as PNG and PDF
ggsave(filename = png_filename, plot = p04, width = 15, height = 10, dpi = 300)  # Save as PNG
ggsave(filename = pdf_filename, plot = p04, width = 15, height = 10)            # Save as PDF


# export cell type abundance data
library(stringr)
library(xlsx)

# export table
tb <- table(colData(sce)$sample_id, colData(sce)$CellType)

write.xlsx(tb,
           file = "results/CellType_abundance.xlsx",
           row.names = TRUE,
           sheetName = str_replace_all("facet by sample id", pattern = "[^[:alnum:]]+", replacement = "_"),
           append = FALSE
)


# Section 5: Find ChAT+ cells markers ---------------------------------------------------------------
library(scran)
library(scater)
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)

# Function for Java and R garbage collection
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}

# Define conditions and cell types
conditions <- levels(colData(sce)$condition)
cell_types <- levels(colData(sce)$CellType)

# Loop over each condition and cell type
for (cond in conditions) {
  # Initialize a list to store DEG data frames for current condition
  DEG_list <- list()
  i <- 0
  for (cell_type in cell_types) {
    i <- i + 1
    # Define clusters of interest (COI)
    COI <- paste(cond, c("Neg", "Pos"), cell_type, sep = "-")
    # Subset SCE object for specific condition and cell type
    sce_subset <- sce[, colData(sce)$txt_ChAT_CellType %in% COI]
    colData(sce_subset)$CellType <- factor(colData(sce_subset)$CellType)
    colData(sce_subset)$ChAT_CellType <- factor(colData(sce_subset)$ChAT_CellType)
    colData(sce_subset)$txt_ChAT_CellType <- factor(colData(sce_subset)$txt_ChAT_CellType)
    levels(colData(sce_subset)$txt_ChAT_CellType)
    # Skip if subset is empty
    if (any(table(colData(sce_subset)$chat_exp) < 100)) next
    # Perform differential expression analysis
    DEG_info <- findMarkers(sce_subset, colData(sce_subset)$txt_ChAT_CellType,
      test.type = "wilcox", 
      assay.type = "counts", BPPARAM = MulticoreParam(nworkers)
    )
    # Save results to list
    cl_name <- paste(cond, "Pos", cell_type, sep = "-")
    cl_DEG <- DEG_info[[cl_name]] %>%
      as.data.frame() %>%
      mutate(AUC_rel = summary.AUC - 0.5) %>%
      mutate(AUC_status = ifelse(AUC_rel > 0, "pos", "neg"))
    file_name <- paste0("results/DEG_", cond, ".xlsx")
    write.xlsx(cl_DEG,
      file = file_name,
      row.names = TRUE,
      sheetName = cl_name,
      append = ifelse(i == 1, FALSE, TRUE)
    )
    rm(cl_DEG)
    jgc()
  }
}

# Section 6: Find Treatment-induced cells DEMs ---------------------------------------------------------------
library(scran)
library(scater)
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)

# Function for Java and R garbage collection
jgc <- function() {
    .jcall("java/lang/System", method = "gc")
    gc()
}

# Define conditions and cell types
conditions <- levels(colData(sce)$condition)
cell_types <- levels(colData(sce)$CellType)

# Loop over each cell type
for (cell_type in cell_types) {
  # Initialize a list to store DEG data frames for current cell type
  DEG_list <- list()
  # Pairwise comparisons
  comparisons <- list(c("LPS", "Con"), c("Pam", "Con"), c("Pam", "LPS"))
  sheet_names <- c("LPS_vs_Con", "Pam_vs_Con", "Pam_vs_LPS")
  # Loop for each comparison
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    # Subset SCE object for specific cell type and conditions
    sce_subset <- sce[, colData(sce)$CellType == cell_type & colData(sce)$condition %in% comp]
    colData(sce_subset)$condition <- factor(colData(sce_subset)$condition)
    colData(sce_subset)$CellType <- factor(colData(sce_subset)$CellType)

    # Skip if subset is empty
    if (any(table(colData(sce_subset)$condition) < 100)) next

    # Perform differential expression analysis
    DEG_info <- findMarkers(sce_subset, colData(sce_subset)$condition,
                            test.type = "wilcox", assay.type = "counts", BPPARAM = MulticoreParam(nworkers))
    
    # Save results to list
    cl_DEG <- DEG_info[[comp[1]]] %>%
      as.data.frame() %>%
      mutate(AUC_rel = summary.AUC - 0.5) %>%
      mutate(AUC_status = ifelse(AUC_rel > 0, "pos", "neg"))
    DEG_list[[sheet_names[i]]] <- cl_DEG

    rm(cl_DEG)
    jgc()
  }

  # Write results to Excel
  file_name <- paste0("results/DEG_", cell_type, ".xlsx")
  for (j in seq_along(DEG_list)) {
    write.xlsx(DEG_list[[j]],
      file = file_name,
      row.names = TRUE,
      sheetName = names(DEG_list)[j],
      append = ifelse(j == 1, FALSE, TRUE)
    )
  }
}

# Section 7: plot Chat+ vs Chat- macrophage DEMs ---------------------------------------------------------------
library(openxlsx)
# read DEM results
DEG_Mf <- read_xlsx("results/DEG_Pam.xlsx", sheet = "Pam-Pos-Mfs")
colnames(DEG_Mf)[1] <- "marker"

View(DEG_Mf)

# plot AUC waterfall 
p05 <- ggplot(DEG_Mf, aes(x = reorder(marker, AUC_rel,decreasing = TRUE), y = AUC_rel, fill = AUC_status)) +
  geom_bar(stat = "identity", position = "identity", colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
  coord_cartesian(ylim = c(-0.5, 0.5), clip = "off") +
    scale_y_continuous(labels = function(x) x + 0.5) +    
  theme(
    axis.ticks.x = element_blank(),             
    axis.text.x = element_blank(),            
    panel.background = element_blank() 
  )
p05
                       
# Define output filenames
png_filename <- "figures/Fig3a_ChAT_Mf_waterfall.png"
pdf_filename <- "figures/Fig3a_ChAT_Mf_waterfall.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p05, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p05, width = 15, height = 10)


# Section 8: Find ChAT+ SPM vs ChAT- SPM DEMs and plot ---------------------------------------------------------------
library(scran)
library(scater)
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)

nworkers <- 16

# Function for Java and R garbage collection
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}

# subset macrophages only
plotDR(sce, "UMAP_type", color_by = "CellType")
names(reducedDims(sce))


UMAP_type_Mf  <- (reducedDims(sce)[["UMAP_type"]][,1] > 3) & (reducedDims(sce)[["UMAP_type"]][,2] < 5)
sce_Mf <- sce[,colData(sce)$CellType %in% c("Mfs") & UMAP_type_Mf]

# see macrophage numbers in each group
table(colData(sce_Mf)$sample_id)

# # Define "type" markers
# rowData(sce_Mf)$marker_class <- "state"
# type_maker_names <- grep("ChAT|Tim_4|CD43|CD24|CD9$|F4_80_APC|SSC_H|CD73|CD31_|CD11a|CCR5|PD_1H", rownames(sce_Mf), value = T)
# type_maker_names

# rowData(sce_Mf)$marker_class[rownames(sce_Mf) %in% type_maker_names] <- "type"
table(rowData(sce_Mf)$marker_class)

# redo dimension reduction, clustering and annotation
library(BiocParallel)
library(scran)
library(bluster)
set.seed(123)

# PCA
sce_Mf <- fixedPCA(sce_Mf, assay.type = "exprs", BPPARAM = MulticoreParam(nworkers)) 
# UMAP
sce_Mf <- runDR(sce_Mf, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworkers)) # features = NULL to use all features
# sce_Mf <- runDR(sce_Mf, "TSNE", features = NULL, BPPARAM = MulticoreParam(nworkers)) # features = NULL to use all features

p06_2 <- plotDR(sce_Mf, "UMAP", color_by = "sample_id")
# Define output filenames
png_filename <- "figures/Fig4a_Mf_samples.png"
pdf_filename <- "figures/Fig4a_Mf_samples.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p06_2, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p06_2, width = 15, height = 10)

# Comparison of markers
p07 <- plot_grid(
  plotDR(sce_Mf, "UMAP", color_by = "ChAT_GFP"),
  plotDR(sce_Mf, "UMAP", color_by = "SSC_H"),
  plotDR(sce_Mf, "UMAP", color_by = "F4_80_APC"),
  plotDR(sce_Mf, "UMAP", color_by = "CD43"),
  plotDR(sce_Mf, "UMAP", color_by = "CD11a"),
  plotDR(sce_Mf, "UMAP", color_by = "CD31_aka_PECAM_1"),
  plotDR(sce_Mf, "UMAP", color_by = "Tim_4"),
  plotDR(sce_Mf, "UMAP", color_by = "CD9"),
  plotDR(sce_Mf, "UMAP", color_by = "CD24")
)

# Define output filenames
png_filename <- "figures/Fig4b_Mf_features.png"
pdf_filename <- "figures/Fig4b_Mf_features.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p07, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p07, width = 15, height = 10)

p08 <- plot_grid(
  plotDR(sce_Mf, "UMAP", color_by = "FSC_H"),
  plotDR(sce_Mf, "UMAP", color_by = "CD73"),
  plotDR(sce_Mf, "UMAP", color_by = "MERTK_aka_Mer"),
  plotDR(sce_Mf, "UMAP", color_by = "CD11b_PE_CF594"),
  plotDR(sce_Mf, "UMAP", color_by = "CD93_aka_AA4_1_early_B_lineage"),
  plotDR(sce_Mf, "UMAP", color_by = "MHC_II_BV421")
)

# Define output filenames
png_filename <- "figures/Fig4C_Mf_features_2.png"
pdf_filename <- "figures/Fig4C_Mf_features_2.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p08, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p08, width = 15, height = 10)


## cluster DEM analysis
# define cluster of interest
UMAP_chat  <- (reducedDims(sce_Mf)[["UMAP"]][,2] > 1) & (colData(sce_Mf)$sample_id == "Pam_GFP_Pos")
head(UMAP_chat)
head(ifelse(UMAP_chat,"Pam3_ChAT", "others"))
colData(sce_Mf)$ChAT_clster <- factor(ifelse(UMAP_chat,"Pam3_ChAT", "others"))
# plot groups
plotDR(sce_Mf, "UMAP", color_by = "ChAT_clster")
# DEG analysis
DEG_info <- findMarkers(sce_Mf, colData(sce_Mf)$ChAT_clster,
      test.type = "wilcox", 
      assay.type = "counts", BPPARAM = MulticoreParam(nworkers)
    )
# Save results
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)
# Function for Java and R garbage collection
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}

cl_DEG <- DEG_info[["Pam3_ChAT"]] %>%
  as.data.frame() %>%
  mutate(AUC_rel = summary.AUC - 0.5) %>%
  mutate(AUC_status = ifelse(AUC_rel > 0, "pos", "neg"))


file_name <- paste0("results/DEG_Pam3_ChAT_Pos_vs_others.xlsx")
write.xlsx(cl_DEG,
  file = file_name,
  row.names = TRUE,
  append = FALSE
)

# watherfall plot
DEG_Mf <- cl_DEG
colnames(DEG_Mf)[1] <- "marker"

View(DEG_Mf)

# plot AUC waterfall 
p08 <- ggplot(DEG_Mf, aes(x = reorder(marker, AUC_rel,decreasing = TRUE), y = AUC_rel, fill = AUC_status)) +
  geom_bar(stat = "identity", position = "identity", colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
  coord_cartesian(ylim = c(-0.5, 0.5), clip = "off") +
    scale_y_continuous(labels = function(x) x + 0.5) +   
  theme(
    axis.ticks.x = element_blank(),              
    axis.text.x = element_blank(),            
    panel.background = element_blank() 
  )
p08
# Define output filenames
png_filename <- "figures/Fig5_ChAT_Mf_vs_others_waterfall.png"
pdf_filename <- "figures/Fig5_ChAT_Mf_vs_others_waterfall.pdf"

# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p08, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p08, width = 15, height = 10)


# ! save data
saveRDS(sce_Mf, "InfinityFlow_Mf-step3_Mf_clustering.rds")
# sce_Mf <- readRDS("InfinityFlow_Mf-step3_Mf_clustering.rds")

# Section 9: Cell cluster DM plots ---------------------------------------------------------------
library(BiocParallel)
library(scran)
library(bluster)
library(rJava)
library(xlsx)
library(stringr)
library(dplyr)
set.seed(123)

# Function for Java and R garbage collection
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
  gc()
}
sce <- readRDS("InfinityFlow_Mf-step3_clustering.rds")
# for subset cell types
plotDR(sce, "UMAP_type", color_by = "CellType")

# macrophages
sce_Mf <- readRDS("InfinityFlow_Mf-step3_Mf_clustering.rds")
# 2-step clustering
primary_cluster  <- min(1000, ncol(sce_Mf))
Mf_clusters <- clusterCells(sce_Mf,
  # assay.type = "exprs", 
  use.dimred="UMAP",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = primary_cluster),
    second = NNGraphParam(k = 40)
  )
)
table(Mf_clusters)
colData(sce_Mf)$Mf_clusters <- Mf_clusters

p_Mf <- plotDR(sce_Mf, "UMAP", color_by = "Mf_clusters")
# Define output filenames
png_filename <- "figures/Mf_clusters.png"
pdf_filename <- "figures/Mf_clusters.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p_Mf, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p_Mf, width = 15, height = 10)
# find feature markers - DEG analysis
DEG_Mf <- scoreMarkers(sce_Mf, colData(sce_Mf)$Mf_clusters,
      assay.type = "exprs", 
      BPPARAM = MulticoreParam(nworkers)
    )
file_path <- "results/DEG_Mf_markers.xlsx"
for (i in seq_along(DEG_Mf)) {
  df <- as.data.frame(DEG_Mf[[i]])
  cluster_name <- names(DEG_Mf)[i]
  sheet_name <- str_replace_all(cluster_name, pattern = "[^[:alnum:]]+", replacement = "_")
  write.xlsx(df, file = file_path, sheetName = sheet_name, append = i != 1, row.names = TRUE)
  jgc()
}
# ! save data
saveRDS(sce_Mf, "InfinityFlow_Mf-step3_Mf_clustering.rds")
# sce_Mf <- readRDS("InfinityFlow_Mf-step3_Mf_clustering.rds")

# B cells
UMAP_type_B  <- (reducedDims(sce)[["UMAP_type"]][,1] < 0) & (reducedDims(sce)[["UMAP_type"]][,2] < 7.5)
sce_B <- sce[,colData(sce)$CellType %in% c("B_cells") & UMAP_type_B]
# PCA
sce_B <- fixedPCA(sce_B, assay.type = "exprs", subset.row=NULL, BPPARAM = MulticoreParam(nworkers)) 
# UMAP
sce_B <- runDR(sce_B, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworkers)) # features = NULL to use all features
# 2-step clustering
primary_cluster  <- min(1000, ncol(sce_B))
B_clusters <- clusterCells(sce_B,
  # assay.type = "exprs", 
  use.dimred="UMAP",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = primary_cluster),
    second = NNGraphParam(k = 50)
  )
)
table(B_clusters)
colData(sce_B)$B_clusters <- B_clusters

p_B <- plotDR(sce_B, "UMAP", color_by = "B_clusters")
# Define output filenames
png_filename <- "figures/B_clusters.png"
pdf_filename <- "figures/B_clusters.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p_B, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p_B, width = 15, height = 10)
# find feature markers - DEG analysis
DEG_B <- scoreMarkers(sce_B, colData(sce_B)$B_clusters,
      assay.type = "exprs", 
      BPPARAM = MulticoreParam(nworkers)
    )
file_path <- "results/DEG_B_markers.xlsx"
for (i in seq_along(DEG_B)) {
  df <- as.data.frame(DEG_B[[i]])
  cluster_name <- names(DEG_B)[i]
  sheet_name <- str_replace_all(cluster_name, pattern = "[^[:alnum:]]+", replacement = "_")
  write.xlsx(df, file = file_path, sheetName = sheet_name, append = i != 1, row.names = TRUE)
  jgc()
}
# ! pause here, save data
saveRDS(sce_B, file = "InfinityFlow_B-step3_clustering.rds")
# sce_B <- readRDS("InfinityFlow_B-step3_clustering.rds")

# T and NK cells
UMAP_type_TNK <- (reducedDims(sce)[["UMAP_type"]][,1] < 5) & (reducedDims(sce)[["UMAP_type"]][,2] > 7.5)
sce_TNK <- sce[,colData(sce)$CellType %in% c("gdT_cells", "abT_cells", "NK_cells") & UMAP_type_TNK]
# PCA
sce_TNK <- fixedPCA(sce_TNK, assay.type = "exprs",subset.row=NULL, BPPARAM = MulticoreParam(nworkers)) 
# UMAP
sce_TNK <- runDR(sce_TNK, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworkers)) # features = NULL to use all features
# 2-step clustering
primary_cluster  <- min(500, ncol(sce_TNK))
TNK_clusters <- clusterCells(sce_TNK,
  # assay.type = "exprs", 
  use.dimred="UMAP",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = primary_cluster),
    second = NNGraphParam(k = 40)
  )
)
table(TNK_clusters)
colData(sce_TNK)$TNK_clusters <- TNK_clusters

p_TNK <- plotDR(sce_TNK, "UMAP", color_by = "TNK_clusters")
# Define output filenames
png_filename <- "figures/T&NK_clusters.png"
pdf_filename <- "figures/T&NK_clusters.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p_TNK, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p_TNK, width = 15, height = 10)
# find feature markers - DEG analysis
DEG_TNK <- scoreMarkers(sce_TNK, colData(sce_TNK)$TNK_clusters,
      assay.type = "exprs", 
      BPPARAM = MulticoreParam(nworkers)
    )
file_path <- "results/DEG_TNK_markers.xlsx"
for (i in seq_along(DEG_TNK)) {
  df <- as.data.frame(DEG_TNK[[i]])
  cluster_name <- names(DEG_TNK)[i]
  sheet_name <- str_replace_all(cluster_name, pattern = "[^[:alnum:]]+", replacement = "_")
  write.xlsx(df, file = file_path, sheetName = sheet_name, append = i != 1, row.names = TRUE)
  jgc()
}
# ! pause here, save data
saveRDS(sce_TNK, file = "InfinityFlow_T_NK-step3_clustering.rds")
# sce_TNK <- readRDS("InfinityFlow_T_NK-step3_clustering.rds")


# myeloid cells (except for Mfs) 
UMAP_type_Mye <- (reducedDims(sce)[["UMAP_type"]][,1] > 3)
sce_Mye <- sce[,colData(sce)$CellType %in% c("DCs", "Neutrophils") & UMAP_type_Mye]
# PCA
sce_Mye <- fixedPCA(sce_Mye, assay.type = "exprs", subset.row=NULL, BPPARAM = MulticoreParam(nworkers)) 
# UMAP
sce_Mye <- runDR(sce_Mye, "UMAP", features = NULL, BPPARAM = MulticoreParam(nworkers)) # features = NULL to use all features
# 2-step clustering
primary_cluster  <- min(500, ncol(sce_Mye))
Mye_clusters <- clusterCells(sce_Mye,
  # assay.type = "exprs", 
  use.dimred="UMAP",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = primary_cluster),
    second = NNGraphParam(k = 40)
  )
)
table(Mye_clusters)
colData(sce_Mye)$Mye_clusters <- Mye_clusters

p_Mye <- plotDR(sce_Mye, "UMAP", color_by = "Mye_clusters")
# Define output filenames
png_filename <- "figures/Mye_clusters.png"
pdf_filename <- "figures/Mye_clusters.pdf"
# Save the combined plot as PNG and PDF
ggsave(filename = png_filename, plot = p_Mye, width = 15, height = 10, dpi = 300)
ggsave(filename = pdf_filename, plot = p_Mye, width = 15, height = 10)
# find feature markers - DEG analysis
DEG_Mye <- scoreMarkers(sce_Mye, colData(sce_Mye)$Mye_clusters,
      assay.type = "exprs", 
      BPPARAM = MulticoreParam(nworkers)
    )
file_path <- "results/DEG_Mye_markers.xlsx"
for (i in seq_along(DEG_Mye)) {
  df <- as.data.frame(DEG_Mye[[i]])
  cluster_name <- names(DEG_Mye)[i]
  sheet_name <- str_replace_all(cluster_name, pattern = "[^[:alnum:]]+", replacement = "_")
  write.xlsx(df, file = file_path, sheetName = sheet_name, append = i != 1, row.names = TRUE)
  jgc()
}

### Final Data Saving --------------------------------------------------------
# Save final processed sce object and analysis results for future use
saveRDS(sce, file = "InfinityFlow_Mf-step3_clustering_20231010.rds")
saveRDS(sce_Mye, file = "InfinityFlow_Mye-step3_clustering.rds")


