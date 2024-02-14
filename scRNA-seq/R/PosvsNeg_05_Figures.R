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
seu <- readRDS("04_WT_Mfs.rds")
seu_SPM <- readRDS("04_WT_SPMs.rds")

cl <- c(
  "LPM_Gata6" = "#EB6F5D",
  "SPM_Itgal" = "#71C9DD",
  "SPM_Lyz1" = "#FA9645",
  "SPM_Cd74" = "#AA8ED6",
  "Mf_Mki67" = "#489746"
)

hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), v = T)
hb_genes <- c(grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), v = T))
bad_features <- unique(c(
    hist_genes, hb_genes,
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
    group.by = "CellType_l2", 
    label = TRUE
) +
    coord_fixed(ratio = 1)
ggsave(file = paste0("figures/WT/fig1_WT_umap_cluster.png"), width = 250, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig1_WT_umap_cluster.pdf"), width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')


# Plot 1: Umap and clusters: WT vs KO
DimPlot(seu,
    reduction = "umap",
    cols = cl,
    split.by = "group",
    group.by = "CellType_l2", 
    ncol = 3,
    label = FALSE
) +
    coord_fixed(ratio = 1)
ggsave(file = paste0("figures/WT/fig1b_WT_split_umap_cluster.png"), width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig1b_WT_split_umap_cluster.pdf"), width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')


### Figure 2: cluster markers --------------------------------------------------
Idents(seu) <- seu$CellType_l2
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

ggsave(file = paste0("figures/WT/fig2a_WT_dotmap_clusterDEGs.png"), width = 100, height = 400, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig2a_WT_dotmap_clusterDEGs.pdf"), width = 100, height = 400, units = "mm", device = "pdf", bg = 'transparent')

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

ggsave(file = paste0("figures/WT/fig2b_WT_dotmap_clusterRegulons.png"), width = 100, height = 400, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig2b_WT_dotmap_clusterRegulons.pdf"), width = 100, height = 400, units = "mm", device = "pdf", bg = 'transparent')


### Figure 3: Cluster abundance in clusters-bar chart --------------------------------------------------
library(tidyverse)
source_cluster <- tibble(seu$CellType_l2,
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
write.csv(source_cluster %>% spread(Cluster, perc), file="results/05_CellType_abundance.csv")

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
ggsave(file = paste0("figures/WT/fig3_cluster_abundance.png"), width = 500, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig3_cluster_abundance.pdf"), width = 500, height = 300, units = "mm", device = "pdf", bg = 'transparent')


### Figure 4: differential genes within SPM clusters: Pos vs Neg  --------------------------------------------------
# load Chat_egfp+ vs Chat_egfp- DE results
load(file = "04_WT_SPM_DEGresults.rds")

# Volcano plot
Volcano_DEG <- ChatDEG_SCT %>%
    tibble::rownames_to_column(var = "symbol")  %>%
    dplyr::filter(!symbol %in% bad_features) %>%
    dplyr::arrange(desc(avg_log2FC))
View(Volcano_DEG)

label_features <- c(
  "Chat", "H2-Eb1", "H2-Ab1", "H2-Aa", "Iglc2", "S100a9", "Cd74", "Stab1", "Neat1", "Mmp9", "Slc13a3", "Folr2", 
  "Cxcl2", "Cxcl13", "Cxcl1", "Pf4", "Ccl2", "Ccrl2", "Ccl4","Slpi","Marco","Alox15","Prg4",
  "Ticam2", "Mfhas1", "Gpr108", "Rps6ka3", "Lrrc14", "Tlr2", "Tirap", "Irf4", "Mapkapk2", "Plcg2", "App", "Otud4", "Tlr1", "Tlr6", "Trim30a", "Irak2"
) %>%  unique()

library(EnhancedVolcano)
ChatDEG_Volcano <- EnhancedVolcano(Volcano_DEG,
                                   lab = Volcano_DEG$symbol,
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
ggsave(file = paste0("figures/WT/fig4a_SPM_DEG_Volcano.png"), width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig4a_SPM_DEG_Volcano.pdf"), width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')


# top regulons
# View(ChatDEG_Bin)
topRegulon <- ChatDEG_Bin %>%
    tibble::rownames_to_column(var = "symbol")  %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    .$symbol

DotPlot(seu_SPM[,seu_SPM$group %in% c("Neg","Pos")],
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

ggsave(file = paste0("figures/WT/fig4b_WT_SPM_DEG_Regulons.png"), width = 100, height = 400, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig4b_WT_SPM_DEG_Regulons.pdf"), width = 100, height = 400, units = "mm", device = "pdf", bg = 'transparent')


### Figure 5: GESA: Pos vs Neg  --------------------------------------------------
# Prepare gene sets
library(fgsea)
library(msigdbr)
set.seed(123)
m_df = msigdbr(species = "Mus musculus")
# Displaying distinct gene set categories and subcategories
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n=23)

# Select GO:BP (Gene Ontology: Biological Process) gene sets
pathwaysDF <- msigdbr("mouse", category="C2", subcategory = "CP:REACTOME")
# Splitting gene symbols by gene set names
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
# Select specific pathways
pathways_select <- c(
  "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES", "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE" , "REACTOME_ERK_MAPK_TARGETS","REACTOME_SIGNALLING_TO_ERKS"
)
# Filtering selected pathways
pathways <- pathways[pathways_select]
# Unlisting genes from selected pathways
Chemical_synapses_genes <- unlist(pathways["REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES"])
TLR2_signaling_genes <- unlist(pathways["REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE"])
ERK_MAPK_targets_genes <- unlist(pathways["REACTOME_ERK_MAPK_TARGETS"])
Signaling_to_ERKS_genes <- unlist(pathways["REACTOME_SIGNALLING_TO_ERKS"])


# Combining all genes into a single data frame
geneset <- data.frame(
  term = c(
    rep("Chemical_synapses", length(Chemical_synapses_genes)),
    rep("TLR2_signaling", length(TLR2_signaling_genes)),
    rep("ERK_MAPK_targets", length(ERK_MAPK_targets_genes)),
    rep("Signaling_to_ERKS", length(Signaling_to_ERKS_genes))
  ),
  gene = c(
    Chemical_synapses_genes,
    TLR2_signaling_genes,
    ERK_MAPK_targets_genes,
    Signaling_to_ERKS_genes
  )
)
# Display the tail of the combined geneset
tail(geneset)

# Prepare gene list for GSEA
load(file = "04_WT_SPM_DEGresults.rds")
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
GESE_res <- GSEA(geneList, TERM2GENE=geneset, verbose=TRUE)

# Plotting GSEA results for the gene sets with adjusted y-axis limits
gseaplot2(GESE_res, geneSetID = 1:4, subplots = 1:2)

ggsave(file = paste0("figures/WT/fig5_WT_GSEA.png"), width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/WT/fig5_WT_GSEA.pdf"), width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')

# Viewing GSEA results
write.csv(GESE_res@result, file="results/05_WT_GSEA.csv")
