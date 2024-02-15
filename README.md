# ChEM Analysis Code Repository
Welcome to the GitHub repository for the analysis code used in the study "Cholinergic Macrophages Promote Resolution of Peritoneal Inflammation," focusing on ChAT-expressing Macrophages (ChEM).

## Overview
This repository contains the computational analysis codes for our research, which elucidates the role of a novel subset of cholinergic macrophages in the resolution of acute peritonitis. By investigating ChAT-expressing macrophages (Mφs), our study highlights the critical role of the non-neural cholinergic system in regulating immune equilibrium and tissue homeostasis.

## Abstract
Our research identifies ChAT-expressing cholinergic macrophages (Mφs) during the resolution phase of acute peritonitis, showcasing their significance in immune regulation. Using Chat-GFP reporter mice, we discovered a marked upregulation of ChAT in monocyte-derived small peritoneal macrophages (SPMs) in response to Pam3CSK4 (Pam3), mediated through a MyD88-dependent MAPK signaling pathway. Functionally, we observed that Chat deficiency in Mφs led to reduced efferocytosis of apoptotic neutrophils and a delayed resolution of peritonitis, reversible with ACh supplementation. These findings underscore the importance of Mφ-derived ACh in inflammation resolution, emphasizing the non-neuronal cholinergic system's role in immune regulation.

## Repository Contents
This repository contains R and Python scripts and analysis pipelines developed and used for the study. The analysis is divided into two main parts:

- **InfinityFlow analysis**: Utilized for analyzing cell surface markers to identify cholinergic macrophages within the peritoneal cavity. This method allows for the high-throughput screening of cell populations based on their phenotypic markers, enabling the distinction of ChEMs from other macrophage subsets.

- **scRNA-seq analysis**: Employed for transcriptomic profiling of the identified cholinergic macrophages, revealing insights into their functional roles during the resolution of inflammation. This comprehensive analysis sheds light on the gene expression patterns that underpin the unique functions of ChEMs, including their contributions to efferocytosis and the modulation of inflammatory responses.

## Getting Started
To use these analysis codes, ensure you have the following R packages installed: Seurat, patchwork, ggplot2, dplyr, future, rJava, xlsx, stringr, SCopeLoomR, AUCell, SCENIC, and other dependencies as mentioned in the scripts.

1. Clone this repository to your local machine.
2. Set your working directory to where the data files are located, as specified in the scripts.
3. Execute the scripts in R or RStudio to replicate the analysis or apply it to your datasets.

# Contact
For any queries related to the code or research, please open an issue on this GitHub repository or email to wuchong5@sysu.edu.cn.
