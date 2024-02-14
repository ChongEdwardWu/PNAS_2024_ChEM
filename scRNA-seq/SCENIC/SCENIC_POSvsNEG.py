import os

# import dependencies
import numpy as np
import pandas as pd
import loompy as lp
from pyscenic.binarization import binarize


# !需要手动设置的变量
DATASET_ID='seu_int'
WORKDIR="/your_work_dir/data/Peritoneal_ChAT"
AUXILLIARIES_FOLDERNAME="/your_work_dir/Ref_genome/SCENIC/Mouse" #! change spesis if necessary!


# set variables for file paths to read from and write to:
R_FOLDERNAME=os.path.join(WORKDIR,"scRNA_seq")
RESULTS_FOLDERNAME=os.path.join(WORKDIR, "SCENIC/results")
FIGURES_FOLDERNAME=os.path.join(WORKDIR, "SCENIC/figures")

# path to loom file converted from the Seurat Object.
INPUTLOOM_FNAME=os.path.join(R_FOLDERNAME,f"{DATASET_ID}.loom")
# Ranking databases. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
RANKING_DBS_FNAMES = " ".join( #将AUXILLIARIES_FOLDERNAME目录下所有结尾为 .rankings.feather 的文件名以空格相连，赋值给 RANKING_DBS_FNAMES 变量。
    [
        os.path.join(AUXILLIARIES_FOLDERNAME, f)
        for f in os.listdir(AUXILLIARIES_FOLDERNAME)
        if f.endswith(".rankings.feather")
    ]
)

# Motif annotations. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
MOTIF_ANNOTATIONS_FNAME=os.path.join(AUXILLIARIES_FOLDERNAME,"motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
MM_TFS_FNAME=os.path.join(AUXILLIARIES_FOLDERNAME,"allTFs_mm.txt")

# Output file
ADJACENCIES_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.adjacencies.tsv")
MOTIFS_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.motifs.csv")
REGULONS_DAT_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.regulons.dat")
AUCELL_MTX_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.auc.csv")
BIN_MTX_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.bin.csv")
THR_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.thresholds.csv")
ANNDATA_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.h5ad")
LOOM_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.scenic.loom")
UMAP_FNAME=os.path.join(RESULTS_FOLDERNAME,f"{DATASET_ID}.umap.csv")

# STEP 1: Gene regulatory network inference, and generation of co-expression modules
if os.system(
    f"pyscenic grn {INPUTLOOM_FNAME} {MM_TFS_FNAME} -o {ADJACENCIES_FNAME} "
    f"--num_workers 16" # ! num_workers
) == 0:
    print("SCENIC STEP 1 is done！")

# STEP 2-3: Regulon prediction aka cisTarget from CLI
if os.system(
    f"pyscenic ctx {ADJACENCIES_FNAME} {RANKING_DBS_FNAMES} "
    f"--annotations_fname {MOTIF_ANNOTATIONS_FNAME} "
    f"--expression_mtx_fname {INPUTLOOM_FNAME} "
    f"--output {MOTIFS_FNAME} "
    f"--auc_threshold 0.05 "
    f"--num_workers 8" # ! num_workers
) == 0:
    print("SCENIC STEP 2 is done！")

if os.system(
    f"pyscenic aucell {INPUTLOOM_FNAME} {MOTIFS_FNAME} "
    f"--output {LOOM_FNAME} " # {AUCELL_MTX_FNAME}
    f"--num_workers 8" # ! num_workers
) == 0:
    print("SCENIC STEP 3 is done！")
    
    
# STEP 4:  Regulon activity binarization
if os.system(
    f"pyscenic aucell {INPUTLOOM_FNAME} {MOTIFS_FNAME} "
    f"--output {AUCELL_MTX_FNAME} " # 
    f"--num_workers 8" # ! num_workers
) == 0:
    print("AUC matrix has been created！")
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)
bin_mtx, thresholds = binarize(auc_mtx,seed=123,num_workers=8) # ! num_workers
bin_mtx.to_csv(BIN_MTX_FNAME) 
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)
print("SCENIC STEP 4 is done！")
