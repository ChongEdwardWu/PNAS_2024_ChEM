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
library(future)
nworkers <- 16
# plan("sequential")
plan("multicore", workers = nworkers)

# Load the processed Seurat object
seu <- readRDS("04_KO_all.rds")
seu_Mf <- readRDS("04_KO_Mfs.rds")

cl <- c(
  "Mf" = "#EB6F5D",
  "B_cell" = "#71C9DD",
  "Neutrophil" = "#FA9645",
  "T_cell" = "#AA8ED6",
  "NK_cell" = "#489746",
  "Mast_cell" = "#8C564B"
)

cl2 <- c(
  "WT" = "#a18ac1",
  "KO" = "#84c9bf"
)

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

### Figure 1: Dimension reduction --------------------------------------------------

# Plot 1: Umap and clusters 
DimPlot(seu,
    reduction = "umap",
    cols = cl,
    #split.by = "group",
    group.by = "CellType_l1", 
    label = TRUE
) +
    coord_fixed(ratio = 1)
ggsave(file = paste0("figures/KO/fig1_KO_umap_cluster.png"), width = 250, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig1_KO_umap_cluster.pdf"), width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')


# Plot 1: Umap and clusters: WT vs KO
DimPlot(seu,
    reduction = "umap",
    cols = cl,
    split.by = "group",
    group.by = "CellType_l1", 
    ncol = 3,
    label = FALSE
) +
    coord_fixed(ratio = 1)
ggsave(file = paste0("figures/KO/fig1b_KO_split_umap_cluster.png"), width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig1b_KO_split_umap_cluster.pdf"), width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')


### Figure 2: cluster markers --------------------------------------------------
Idents(seu) <- seu$CellType_l1
topDEG <- unique(c(
  "Gata6","Prg4", "Itga6", "Cd9", "Pf4", "Arg1", "Marco", "Cd24a", 
  "Itgal", "Efhd2", "Kctd12", "Tmsb4x", "Lst1", "Lcp1", "Anxa5", "Ptpn1", "Aif1", "Pfn1", "Pecam1", 
  "Lyz1", "Alox5ap", "Cd300c2", "Ifitm2", "H3f3a", 
  "Cd74", "H2-Aa", "H2-Eb1", "H2-Ab1", "Clec4b1", 
  "Mki67", "Birc5", "Stmn1", "Ccnb1", "Cenpa"
))

DotPlot(seu,
    assay = "SCT",
    # group.colors = cl,
    c("white", "#AD272B"),
    features = rev(topDEG)
) +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "bottom")

ggsave(file = paste0("figures/KO/fig2a_KO_dotmap_clusterDEGs.png"), width = 100, height = 400, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig2a_KO_dotmap_clusterDEGs.pdf"), width = 100, height = 400, units = "mm", device = "pdf", bg = 'transparent')

topRegulon <- unique(c(
  "Gata6", "Nfia",  "Xbp1", "Hes1", 
  "Spi1", "Nfatc2", "Klf4", "Arid3a", 
  "Ehf", "Hoxa7", "Runx3", "Sox4", 
  "Gabpa", "E2f8", "Brca1", "Nfyb", "Sp2"
))

DotPlot(seu,
    assay = "Bin",
    # group.colors = cl,
    c("white", "#AD272B"),
    features = rev(topRegulon)
) +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "bottom")

ggsave(file = paste0("figures/KO/fig2b_KO_dotmap_clusterRegulons.png"), width = 100, height = 400, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig2b_KO_dotmap_clusterRegulons.pdf"), width = 100, height = 400, units = "mm", device = "pdf", bg = 'transparent')


### Figure 3: Cluster abundance in clusters-bar chart --------------------------------------------------
library(tidyverse)
source_cluster <- tibble(seu$CellType_l1,
                              seu$group) %>%
  set_names("Cluster", "Group") %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n()) %>%
  group_by(Group) %>%
  mutate(total.no = sum(no.cell),
         perc = 100*no.cell/total.no) %>%
  dplyr::select(Cluster, Group, perc) 
View(source_cluster %>%
  spread(Cluster, perc))
write.csv(source_cluster %>% spread(Cluster, perc), file="results/05_KO_CellType_abundance.csv")

library(ggplot2)
View(source_cluster)
ggplot(
  source_cluster,
  aes(x = Group, y = perc, fill = Cluster)
) +
  geom_col(colour = "black") +
  scale_fill_manual(values = cl) +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("") +
  ylab("%")
ggsave(file = paste0("figures/KO/fig3_KO_cluster_abundance.png"), width = 500, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig3_KO_cluster_abundance.pdf"), width = 500, height = 300, units = "mm", device = "pdf", bg = 'transparent')


### Figure 4: differential genes within clusters (except for macrophages): KO vs WT  --------------------------------------------------
# Load necessary libraries
library(Seurat)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(tibble)
library(xlsx)

# Function to plot and save volcano plots
plot_save_volcano <- function(deg_data, cell_type) {
    # Prepare data for volcano plot
    volcano_data <- deg_data %>%
        rownames_to_column(var = "symbol") %>%
        filter(!symbol %in% bad_features) %>%
        arrange(p_val_adj)
    
    # Select top 20 genes
    selected_genes <- volcano_data$symbol[1:20]

    # Create volcano plot
    volcano_plot <- EnhancedVolcano(
      volcano_data,
      title = paste0(cell_type),
      lab = volcano_data$symbol,
      x = "avg_log2FC",
      y = "p_val_adj",
      xlim = c(-3, 3),
      selectLab = selected_genes,
      pCutoff = 0.05,
      FCcutoff = 0.25,
      xlab = bquote(~ Log[2] ~ "fold change"),
      pointSize = 1.0,
      labSize = 6.0,
      labCol = "black",
      labFace = "bold",
      boxedLabels = TRUE,
      colAlpha = 4 / 5,
      legendPosition = "right",
      legendLabSize = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      lengthConnectors = 2,
      arrowheads = FALSE,
      colConnectors = "black"
    )

    # Save the plot as PNG and PDF
    png_filename <- paste0("figures/KO/fig4_KO_", cell_type, "_DEG_Volcano.png")
    pdf_filename <- paste0("figures/KO/fig4_KO_", cell_type, "_DEG_Volcano.pdf")
    
    ggsave(file = png_filename, plot = volcano_plot, width = 300, height = 300, units = "mm", dpi = 300, device = "png")
    ggsave(file = pdf_filename, plot = volcano_plot, width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')
}

# Define cell types to analyze
cell_types <- setdiff(levels(seu$CellType_l1),"Mf")

# Main loop for each cell type
for (cell_type in cell_types) {
    # Load differential expression results
    result_file_name <- paste0("04_KO_", cell_type, "_DEGresults_", Sys.Date(), ".rds")
    load(file = result_file_name)

    # Call the plotting function
    plot_save_volcano(ChatDEG_SCT, cell_type)
    # Clear variables before next iteration (optional, for memory management)
    rm(ChatDEG_SCT, ChatDEG_AUC, ChatDEG_Bin)
}



### Figure 5: differential genes within macrophages clusters: KO vs WT  --------------------------------------------------
# load Chat_WT+ vs Chat_KO DE results
load(file = "04_KO_Mf_DEGresults.rds")

# Volcano plot
# Prepare data for volcano plot
volcano_data <- ChatDEG_SCT %>%
  rownames_to_column(var = "symbol") %>%
  filter(!symbol %in% bad_features) %>%
  arrange(p_val_adj)

label_features <- c(
  "Alox15", "Tgfb2", "Ltbp1", "Cxcl13",  "Aspa", "Prg4", "Garnl3", "Timd4", "F5", "Plxdc2", 
  "Lyz1", "Lyz2", 
  "Chat", "Irf1", "Irgm2","Stat1",
  volcano_data$symbol[1:20]
) %>%  unique()

library(EnhancedVolcano)
ChatDEG_Volcano <- EnhancedVolcano(volcano_data,
                                   lab = volcano_data$symbol,
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   xlim = c(-3, 3),
                                   selectLab = label_features ,
                                   pCutoff = 0.05,
                                   FCcutoff = 0.25,
                                   xlab = bquote(~Log[2]~ 'fold change'),
                                   pointSize = 1.0,
                                   labSize = 6.0,
                                   labCol = 'black',
                                   labFace = 'bold',
                                   boxedLabels = T,
                                   colAlpha = 4/5,
                                   legendPosition = 'right',
                                   legendLabSize = 14,
                                   legendIconSize = 4.0,
                                   drawConnectors = T,
                                   widthConnectors = 0.5,
                                   lengthConnectors = 2,
                                   arrowheads = F,
                                   colConnectors = 'black')
ChatDEG_Volcano
ggsave(file = paste0("figures/KO/fig5_KO_Mf_DEG_Volcano.png"), width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig5_KO_Mf_DEG_Volcano.pdf"), width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')



### Figure 6: GESA: KO vs WT  --------------------------------------------------
# Prepare gene sets
library(fgsea)
library(msigdbr)
set.seed(123)
m_df = msigdbr(species = "Mus musculus")
# Displaying distinct gene set categories and subcategories
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n=23)

# Select GO:BP (Gene Ontology: Biological Process) gene sets
pathwaysDF <- msigdbr("mouse", category="C5", subcategory = "GO:BP")
# Splitting gene symbols by gene set names
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
# Select specific pathways
pathways_select <- c(
  "GOBP_REGULATION_OF_ENDOCYTOSIS", "GOBP_APOPTOTIC_CELL_CLEARANCE","GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE"
)
# Filtering selected pathways
pathways <- pathways[pathways_select]
# Unlisting genes from selected pathways
Endocytosis_genes <- unlist(pathways["GOBP_REGULATION_OF_ENDOCYTOSIS"])
Apoptotic_cell_clearance_genes <- unlist(pathways["GOBP_APOPTOTIC_CELL_CLEARANCE"])
Innate_immune_response_genes <- unlist(pathways["GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE"])

# Resolution phase Macrophage gene set
Resolution_phase_Macrophage_genes <- c(
  "Alox15", "Gpb1", "Tgfb2", "Ltbp1", "Cxcl13", "Selp", "Aspa", "Prg4", "Garnl3", "Timd4", "F5", "St8sia6",
  "Plxdc2", "Clec1a", "Epha4", "Sele", "Sult1a1", "Ube2c", "Cnnb2", "Shcbp1"
) 

# Combining all genes into a single data frame
geneset <- data.frame(
  term = c(
    rep("Endocytosis", length(Endocytosis_genes)),
    rep("Apoptotic_cell_clearance", length(Apoptotic_cell_clearance_genes)),
    rep("Innate_immune_response", length(Innate_immune_response_genes)),
    rep("Resolution_phase_Macrophage", length(Resolution_phase_Macrophage_genes))
  ),
  gene = c(
    Endocytosis_genes,
    Apoptotic_cell_clearance_genes,
    Innate_immune_response_genes,
    Resolution_phase_Macrophage_genes
  )
)
# Display the tail of the combined geneset
tail(geneset)

# Prepare gene list for GSEA
load(file = "04_KO_Mf_DEGresults.rds")
# Extracting log fold change values and assigning gene names
logFC <- ChatDEG_SCT$avg_log2FC
names(logFC) <- rownames(ChatDEG_SCT)
# Sorting gene list in decreasing order of logFC
geneList <- sort(logFC, decreasing = T)
# Removing genes without names
geneList <- geneList[!is.na(names(geneList))]
# Display the length of the gene list
length(geneList)

# Perform GSEA
library(clusterProfiler)
library(fgsea)
library(enrichplot)
# Running GSEA using the prepared gene list and gene set
rM_GSEA <- GSEA(geneList, TERM2GENE=geneset, verbose=TRUE)

# Plotting GSEA results for the gene sets with adjusted y-axis limits
gseaplot2(rM_GSEA, geneSetID = 1:4, subplots = 1:2)

ggsave(file = paste0("figures/KO/fig6_KO_GSEA.png"), width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig6_KO_GSEA.pdf"), width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')

# Viewing GSEA results
write.csv(rM_GSEA@result, file="results/05_KO_GSEA_results.csv")


### Figure 7: Cell Cycle: KO vs WT  --------------------------------------------------
# cell cycle phages data
cell_cycle_data <- as.matrix(table(seu_Mf$group,seu_Mf$CCphase))

# Perform Chi-square test
chisq_test_result <- chisq.test(cell_cycle_data)
print(chisq_test_result)

# Convert matrix to data frame for easier manipulation
cell_cycle_df <- as.data.frame(cell_cycle_data)

# Rename columns and transform data
cell_cycle_tidy <- cell_cycle_df %>%
  set_names("Group", "Phase", "Count") %>%
  group_by(Group) %>%
  mutate(TotalCount = sum(Count),
         Percent = 100 * Count / TotalCount) %>%
  dplyr::select(Group, Phase, Percent)
write.csv(cell_cycle_tidy , file="results/05_KO_Mf_CellCycle.csv")

# Plotting
ggplot(cell_cycle_tidy, aes(x = Group, y = Percent, fill = rev(Phase))) +
  geom_col(colour = "black") +
  scale_fill_brewer(palette = "Set1") + # You can change the color palette
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("") +
  ylab("% Cell cycle phase")

# save the plot
ggsave(file = paste0("figures/KO/fig7_KO_Mf_Cellcycle.png"), width = 150, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig7_KO_Mf_Cellcycle.pdf"), width = 150, height = 300, units = "mm", device = "pdf", bg = 'transparent')

### Figure 8: RNA splicing: KO vs WT  --------------------------------------------------

# Extracting spliced and unspliced RNA data
spliced_data <- LayerData(seu_Mf, assay = "spliced", layer = "data")
unspliced_data <- LayerData(seu_Mf, assay = "unspliced", layer = "data")

# Calculating total spliced and unspliced RNA for each cell
spliced_totals <- Matrix::colSums(spliced_data)
unspliced_totals <- Matrix::colSums(unspliced_data)

# Calculating the percentage of spliced and unspliced RNA
total_expression <- spliced_totals + unspliced_totals
spliced_percents <- spliced_totals / total_expression * 100
unspliced_percents <- unspliced_totals / total_expression * 100

# Combining data
rna_data <- data.frame(
  Group = seu_Mf$group,
  SplicedPercent = spliced_percents,
  UnsplicedPercent = unspliced_percents
)

# Performing statistical tests
wilcox_result <- wilcox.test(UnsplicedPercent ~ Group, data = rna_data)
p_value <- wilcox_result$p.value
p_label <- paste("Wilcox p-value:", format(p_value, digits = 2))

# Violin plot for UnsplicedPercent
ggplot(rna_data, aes(x = Group, y = UnsplicedPercent, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", 
               outlier.shape = 21, 
               outlier.colour = "black", 
               outlier.fill = "white", 
               outlier.size = 5) +
  scale_fill_manual(values = cl2) +
  labs(title = "Unspliced RNA Percentage by Group",
       x = "Group",
       y = "Unspliced RNA Percentage (%)") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(rna_data$UnsplicedPercent) * 0.9, label = p_label, size = 4)

# save the plot
ggsave(file = paste0("figures/KO/fig8_KO_Mf_UnsplicedPercent.png"), width = 150, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig8_KO_Mf_UnsplicedPercent.pdf"), width = 150, height = 300, units = "mm", device = "pdf", bg = 'transparent')

### Figure 9: Chat expression: KO vs WT  --------------------------------------------------
seu$group_celltype <- factor(seu$group_celltype, levels = c("Mf_WT", "Mf_KO", "B_cell_WT", "B_cell_KO", 
                                                            "T_cell_WT", "T_cell_KO", "NK_cell_WT", 
                                                            "NK_cell_KO", "Neutrophil_WT", "Neutrophil_KO", 
                                                            "Mast_cell_WT", "Mast_cell_KO"))
levels(seu$group_celltype)
DotPlot(seu,
    assay = "SCT",
    features = "Chat",
    group.by = "group_celltype", 
    c("white", "#AD272B")
) +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "right")
ggsave(file = paste0("figures/KO/fig9_KO_Mf_dotmap_Chat.png"), width = 150, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/KO/fig9_KO_Mf_dotmap_Chat.pdf"), width = 150, height = 150, units = "mm", device = "pdf", bg = 'transparent')
