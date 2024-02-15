### Section 0: Preparation --------------------------------------------------

# clear the environment
rm(list = ls())
graphics.off()
gc()

### Section 1-Sepcify sample name and input ------------------------------------------

## ## ## ## ## ## ## ## ## ## ## ##
##  Reading in 10x Genomics data ##
## ## ## ## ## ## ## ## ## ## ## ##
set.seed(123)
nworkers <- 32

## !input sample name
sample <- "Pos"

# set your work diretory
workdir <- file.path("/ifs1/User/yangbinVM/01.Project/wuchong/data/scRNAseq-R", sample)
if (file.exists(workdir)) {
} else {
  dir.create(workdir)
}
if (file.exists(file.path(workdir,"figures"))) {
} else {
  dir.create(file.path(workdir,"figures"))
}
setwd(workdir)

# Redirect output and messages to the log file
current_date <- format(Sys.Date(), "%Y%m%d")
log_file <- paste0(workdir,"/",sample,"_scRNA_QC_step1_", current_date, ".log")
output_connection <- file(log_file, open = "wt")
sink(output_connection)
sink(output_connection, type = "message")

## !specify the CRcount result here
suppressMessages(suppressWarnings(library(DropletUtils)))
suppressMessages(suppressWarnings(library(BiocParallel)))
sce <- read10xCounts("/ifs1/User/yangbinVM/01.Project/wuchong/data/Peritoneal_ChAT/CRcount/Pos/outs/filtered_feature_bc_matrix", col.names = T)

### Section 2-Preprocess before QC ------------------------------------------

## ## ## ## ## ## ## ## ## ## ## ## ##
##  Cell calling for droplet data   ##
## ## ## ## ## ## ## ## ## ## ## ## ##
a <- colSums(counts(sce))
b <- sort(a)[2] # find the second least UMI count as the lower bound
# emptyDrops
e.out <- emptyDrops(counts(sce), lower = b, BPPARAM = MulticoreParam(nworkers))
## Check cell calling for droplet data
summary(e.out$FDR <= 0.001)

# plot mito vs detected features
png(file = "figures/01_Cell_Calling_for_Droplets.png", width = 250, height = 250, units = "mm", res = 150)
plot(e.out$Total, -e.out$LogProb,
  col = ifelse(e.out$FDR <= 0.001, "red", "black"),
  xlab = "Total UMI count", ylab = "-Log Probability"
)
dev.off()

# The Limited field in the output indicates whether or not the computed
# p-value for a particular barcode is bounded by the number of iterations.
# If any non-significant barcodes are TRUE for Limited,
# we may need to increase the number of iterations.
table(Sig = e.out$FDR <= 0.001, Limited = e.out$Limited)

# subset our SingleCellExperiment object to retain only the detected cells
sce <- sce[,which(e.out$FDR <= 0.001)]

## ## ## ## ## ## ## ## ## ## ## ## ##
##  Normalization by deconvolution  ##
## ## ## ## ## ## ## ## ## ## ## ## ##

# We use a pre-clustering step with quickCluster()
# where cells in each cluster are normalized separately
# and the size factors are rescaled to be comparable across clusters.
suppressMessages(suppressWarnings(library(scran)))
set.seed(123)

clust <- quickCluster(sce, BPPARAM = MulticoreParam(nworkers))
table(clust)
deconv.sf <- calculateSumFactors(sce,
  cluster = clust,
  BPPARAM = MulticoreParam(nworkers)
)
summary(deconv.sf)

# Applying the size factors
# Scaling and log-transforming
suppressMessages(suppressWarnings(library(scater)))
norm <- computeSumFactors(sce, cluster = clust, min.mean = 0.1, BPPARAM = MulticoreParam(nworkers))
norm <- logNormCounts(norm)

## ## ## ## ## ## ## ## ##
##   Feature selection  ##
## ## ## ## ## ## ## ## ##

# Compute the variance of the log-normalized expression values
dec <- modelGeneVar(norm)

# Visualizing the fit:
fit <- metadata(dec)

png(file = "figures/02_Feature_selection.png",width = 250, height = 250, units = "mm", res = 150) 
par(mfrow = c(1, 1))
plot(fit$mean, fit$var,
  xlab = "Mean of log-expression",
  ylab = "Variance of log-expression", main = ""
)
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
dev.off()

# Selecting highly variable genes
hvg.var <- getTopHVGs(dec, n = 2000)
length(hvg.var)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##   Dimensionality Reduction and Clustering ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# PCA
suppressMessages(suppressWarnings(library(BiocSingular)))
set.seed(123)
norm <- fixedPCA(norm,
  subset.row = hvg.var,
  BSPARAM = RandomParam(), name = "PCA"
)

# UMAP
norm <- runUMAP(norm,
  dimred = "PCA",
  BPPARAM = MulticoreParam(nworkers)
)

# clustering
suppressMessages(suppressWarnings(library(scran)))
suppressMessages(suppressWarnings(library(bluster)))

colLabels(norm) <- clusterCells(norm,
  use.dimred = "PCA",
  BLUSPARAM = NNGraphParam(type = "jaccard")
)

# plot the clusters
suppressMessages(suppressWarnings(library(ggplot2)))
sceumap <- plotReducedDim(norm, "UMAP", colour_by = "label", text_by = "label")
ggsave(
  filename = "figures/03_Dimensionality_Reduction_and_Clustering.png",
  plot = sceumap,
  width = 250, height = 250, units = "mm",
  dpi = 150, device = "png", bg = "white"
)

## ## ## ## ## ## ## ## ## ## ## ## ##
## Doublet detection by simulation  ##
## ## ## ## ## ## ## ## ## ## ## ## ##
suppressMessages(suppressWarnings(library(scDblFinder)))
suppressMessages(suppressWarnings(library(scran)))
set.seed(123)

# Setting up the parameters for consistency with denoisePCA();
# this can be changed depending on your feature selection scheme.
dbl <- scDblFinder(
  norm,
  clusters = clust,
  nfeatures = 2000,
  BPPARAM = MulticoreParam(nworkers)
)
table(dbl$scDblFinder.class)
norm$DoubletScore <- dbl$scDblFinder.score
norm$DoubletClass <- dbl$scDblFinder.class

# Section 3: Preliminary QC -----------------------------------------------
suppressMessages(suppressWarnings(library(scater)))
set.seed(123)

## ## ## ## ## ## ## ## ## ## ## 
##  Quality Control - Part 2  ##
#     Mark abnormal 'cells'   ##
## ## ## ## ## ## ## ## ## ## ## 

# Retrieving the mitochondrial transcripts using genomic locations included in
# the row-level annotation for the SingleCellExperiment.
is.mito <- grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", rowData(norm)$Symbol)
is.rp <- grep("^Rp[sl]|^RP[SL]", rowData(norm)$Symbol)
is.hsp <- grep("^HSP|^DNAJ|^Hsp|^Dnaj", rowData(norm)$Symbol)
is.hemoglobin <- grep("^Hb[ab]-|^HB[^(P)]", rowData(norm)$Symbol)
is.plat <- grep("Pecam1|Pf4|PECAM1|PF4", rowData(norm)$Symbol)

qc <- perCellQCMetrics(norm, subsets = list(
  Mito = is.mito, Heatshock = is.hsp, Rp = is.rp,
  RBC = is.hemoglobin, Plat = is.plat
))
colData(norm) <- cbind(colData(norm), qc)
colnames(colData(norm))

# (1) Filtering with fixed thresholds
qc.nexprs <- qc$detected < 300
qc.mito <- qc$subsets_Mito_percent > 10
qc.rbc <- qc$subsets_RBC_percent > 1
fix.discard <- qc.nexprs | qc.mito | qc.rbc
table(fix.discard)
norm$fix.discard <- fix.discard

# # (2) Alternatively,
# # We can use methods from robustbase to quantify the "outlyingness" of each
# # cells based on their QC metrics, and then use isOutlier() to identify
# # low-quality cells that exhibit unusually high levels of outlyingness.
# stats <- cbind(
#   log10(qc$sum),
#   log10(qc$detected),
#   qc$subsets_Mito_percent
# )
# suppressMessages(suppressWarnings(library(robustbase)))
# outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
# multi.outlier <- isOutlier(outlying, type = "higher", nmads = 5)
# # attr(multi.outlier, "thresholds")
# table(multi.outlier)
# norm$adapt.discard <- multi.outlier

# table(cbind(fix.discard, multi.outlier))

# plot mito vs detected features
png(file = "figures/04_Cell_Filtering.png",width = 250, height = 250, units = "mm", res = 150) 
plot(qc$detected, qc$subsets_Mito_percent,
  log = "x",
  xlab = "Detected features", ylab = "Mitochondrial %"
)
points(qc$detected[fix.discard], qc$subsets_Mito_percent[fix.discard],
  col = "dodgerblue", pch = 16, cex = 0.5
)
# points(qc$sum[multi.outlier], qc$subsets_Mito_percent[multi.outlier],
#   col = "red", pch = 16, cex = 0.5
# )
legend("topright",
  legend = c("Fix", "Adaptive"),
  col = c("dodgerblue", "red"), lty = 2, cex = 1
)
dev.off()

discard <- fix.discard
# discard <- fix.discard | multi.outlier
# summary(cbind(fix.discard, multi.outlier, discard))
norm$discard <- discard

# Check diagnostic plots
diagnostic <- gridExtra::grid.arrange(
  plotColData(norm, x = "DoubletClass", y = "sum", colour_by = "DoubletClass") +
    ggtitle("Doublets"),
  plotColData(norm, x = "discard", y = "sum", colour_by = "discard") +
    ggtitle("Total count"),
  plotColData(norm,
    x = "discard", y = "subsets_Mito_percent",
    colour_by = "discard"
  ) + ggtitle("Mito percent"),
  plotColData(norm,
    x = "discard", y = "subsets_Rp_percent",
    colour_by = "discard"
  ) + ggtitle("ribosomal protein percent"),
  plotColData(norm,
    x = "discard", y = "subsets_Heatshock_percent",
    colour_by = "discard"
  ) + ggtitle("Heatshock protein percent"),
  plotColData(norm, x = "sum", y = "detected", colour_by = "discard") +
    ggtitle("Detected features"),
    plotColData(norm,
    x = "sum", y = "subsets_RBC_percent",
    colour_by = "discard"
  ) + ggtitle("hemoglobin percent"),
  plotColData(norm,
    x = "sum", y = "subsets_Plat_detected",
    colour_by = "discard"
  ) + ggtitle("platelet percent"),
  ncol = 3
)
ggsave(
  filename = "figures/05_Diagnostic_Plots.png",
  plot = diagnostic,
  width = 500, height = 250, units = "mm",
  dpi = 150, device = "png", bg = "white"
)

# Check cluster quality
clusQC <- gridExtra::grid.arrange(
  plotColData(norm, x = "label", y = "DoubletScore", colour_by = "DoubletClass") +
    ggtitle("Doublets"),
  plotColData(norm, x = "label", y = "sum", colour_by = "label") +
    ggtitle("Total count"),
  plotColData(norm, x = "label", y = "detected", colour_by = "label") +
    ggtitle("Detected features"),
  plotColData(norm, x = "label", y = "subsets_Mito_percent", colour_by = "label") +
    ggtitle("Mito percent"),
  plotColData(norm, x = "label", y = "subsets_Rp_percent", colour_by = "label") +
    ggtitle("ribosomal genes percent"),
  plotColData(norm, x = "label", y = "subsets_Heatshock_percent", colour_by = "label") +
    ggtitle("heat shock protein genes percent"),
  plotColData(norm, x = "label", y = "subsets_RBC_percent", colour_by = "label") +
    ggtitle("RBC percent"),
  plotColData(norm, x = "label", y = "subsets_Plat_percent", colour_by = "label") +
    ggtitle("platelets percent"),
  plotColData(norm, x = "label", y = "discard", colour_by = "label") +
    ggtitle("Discard"),
  ncol = 3
)
ggsave(
  filename = "figures/06_Cluster_QC_Plots.png",
  plot = clusQC,
  width = 500, height = 500, units = "mm",
  dpi = 150, device = "png", bg = "white"
)

# see preliminary filtering results in clusters
sceumap <- gridExtra::grid.arrange(
  plotReducedDim(norm[, norm$discard == 0], "UMAP", colour_by = "label", text_by = "label") +
    ggtitle("Cells remained"),
  plotReducedDim(norm[, norm$discard == 1], "UMAP", colour_by = "label", text_by = "label") +
    ggtitle("Cells discarded"),
  ncol = 2
)
ggsave(
  filename = "figures/07_PreFiltering_UMAP.png",
  plot = sceumap,
  width = 250, height = 250, units = "mm",
  dpi = 150, device = "png", bg = "white"
)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# pause here
setwd(workdir)
saveRDS(norm, file = paste0("01_QC_step1_",sample,".rds"))
print(paste0("Job ",sample," has been sucessfully done at ", Sys.time(), "."))
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
sink()
sink(type = "message")
rm(list = ls())
graphics.off()
gc()

# norm <- norm[, norm$DoubletClass == "singlet"]

# !Run in server
# nohup Rscript --vanilla /ifs1/User/yangbinVM/01.Project/wuchong/data/scRNAseq-R/Pos/Pos_scRNA_QC_step1.R &
