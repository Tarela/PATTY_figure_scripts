# PATTY Manuscript Source Data and Figure Generation Code

This repository contains the source data and analysis/plotting scripts used in the PATTY manuscript.

PATTY is a computational framework designed to correct open-chromatin bias in CUT&Tag data, particularly for histone modifications such as H3K27me3 and H3K27ac. The repository mainly focuses on reproducing manuscript figures, supplementary figures, and associated source data files.

---

# Repository Structure

```text
.
├── F3_crossValidation/
├── F7_scCnT/
├── F146_promoterPattern/
├── SF5_bivalent/
├── SF11_otherHM/
├── SF456_trueFalseRegionPattern/
└── sourceData/
```

---

# Folder Description

## `F3_crossValidation/`

Scripts and source data for cross-validation benchmarking across machine learning and deep learning models.

Includes:
- Logistic Regression (LR)
- Random Forest (RF)
- Gradient Boosting Machine (GBM)
- CNN
- MLP
- RNN
- GRU

Main analyses:
- Model accuracy comparison
- Correlation with gene expression
- Correlation with histone modification profiles
- Input feature benchmarking
- Boxplot visualization of cross-validation performance

Related manuscript figures:
- Figure 3
- Supplementary Figure 2

---

## `F7_scCnT/`

Analysis of single-cell CUT&Tag (nanoCT) datasets.

Includes:
- UMAP visualization
- Cluster comparison
- Confusion matrices
- PATTY-corrected vs raw signal clustering
- H3K27me3 and H3K27ac single-cell analyses

Related manuscript figures:
- Figure 7

---

## `F146_promoterPattern/`

Promoter-centric signal pattern analyses around transcription start sites (TSSs).

Includes:
- Aggregate signal profiles
- Heatmaps
- CUT&Tag vs ChIP-seq comparisons
- Active vs repressed gene analyses
- RNAPII/GRO-seq comparisons
- High-salt CUT&Tag analyses

Related manuscript figures:
- Figures 1, 4, and 6

---

## `SF5_bivalent/`

Analysis of bivalent chromatin regions.

Includes:
- Venn diagram generation
- Overlap comparison between PATTY and MACS2 peak calls
- Bivalent domain analyses

Related manuscript figures:
- Supplementary Figure 5

---

## `SF11_otherHM/`

Analysis of additional histone modifications using pretrained PATTY models.

Includes:
- H3K4me1
- H3K4me3
- H3K36me3

Main analyses:
- Correlation with nearby gene expression
- Comparison between raw CUT&Tag, IgG-normalized, and PATTY-corrected signals

Related manuscript figures:
- Supplementary Figure 11

---

## `SF456_trueFalseRegionPattern/`

Analysis of signal patterns in true-positive and false-positive genomic regions.

Includes:
- Aggregate TSS profiles
- Heatmaps
- CUT&Tag vs ChIP-seq comparisons
- H3K27me3 and H3K27ac analyses

Related manuscript figures:
- Supplementary Figures 4–6

---

## `sourceData/`

Contains source data files corresponding to manuscript figures and supplementary figures.

Example files:
- `F3_b_LR_sourceData.txt`
- `F3_c_CNN_sourceData.txt`
- `F7_a-c_sourceData.txt`
- `F4_A-F_sourceData.txt`

These files were used to generate the figures shown in the manuscript.

---

# Software Requirements

Main analyses were performed in R.

Required R packages include:

```r
ggplot2
dplyr
tidyr
tidyverse
Seurat
mclust
cluster
gridExtra
VennDiagram
ggpattern
```

Some analyses additionally require:
- BEDTools
- standard genomics utilities
- preprocessed CUT&Tag / ChIP-seq / ATAC-seq signal matrices

---

# Data Availability

All source data used to generate manuscript figures are provided in the `sourceData/` directory or generated directly by the included scripts.

The repository contains processed matrices and figure-level source data used for:
- cross-validation benchmarking,
- promoter/TSS signal profiling,
- single-cell CUT&Tag clustering analyses,
- bivalent chromatin analyses,
- and histone modification benchmarking.

Raw sequencing datasets were obtained from publicly available GEO datasets as described in the manuscript.

---

# Reproducibility

Most figures can be reproduced directly from the provided scripts after adjusting local file paths.

Typical workflow:
1. Prepare processed signal matrices or source data files
2. Run figure-specific R scripts
3. Export figures as PDF/PNG

Many scripts directly generate publication-quality figures used in the manuscript.

---

# Figure Naming Convention

Repository folder names follow the manuscript figure structure:

| Folder | Corresponding Figures |
|---|---|
| `F3_crossValidation` | Figure 3 |
| `F7_scCnT` | Figure 7 |
| `F146_promoterPattern` | Figures 1, 4, 6 |
| `SF5_bivalent` | Supplementary Figure 5 |
| `SF11_otherHM` | Supplementary Figure 11 |
| `SF456_trueFalseRegionPattern` | Supplementary Figures 4–6 |

---

# Example Usage

Example R execution:

```r
source("F3_crossValidation/F3_crossValidation.R")
```

Example output:
- `F3B_raw.pdf`
- `F3C_raw.pdf`
- `F3D_withdots.pdf`
- `F3E_withdots.pdf`

---

# Color Scheme

Several figures use the following PATTY-associated color palette:

```r
colors <- c(
  "#DD864E",
  "#969696",
  "#E6E6E6",
  "#C6DBEF",
  "#6BAED6",
  "#2171B5",
  "#08306B"
)
```

The orange color (`#DD864E`) is consistently used to represent PATTY-corrected signals throughout the manuscript.

---

# Single-cell Analyses

Single-cell CUT&Tag analyses were performed using:
- Seurat
- UMAP embeddings
- clustering consistency analyses
- confusion matrix evaluation

The repository includes:
- processed clustering labels,
- UMAP coordinates,
- and plotting scripts for reproducing the manuscript figures.

---

# Genomic Signal Visualization

Promoter and TSS-centered analyses include:
- aggregate signal profiles,
- heatmaps,
- CUT&Tag vs ChIP-seq comparisons,
- active/repressed gene analyses,
- and accessibility-associated false-positive pattern analyses.

Most signal visualization functions are implemented directly in base R.

---

# Notes

- Most scripts are figure-oriented and were used directly for manuscript figure generation.
- File paths in scripts may need adjustment depending on local directory structure.
- Some large intermediate files and raw sequencing datasets are not included in this repository.

---

# Disclaimer

This repository is intended for research and reproducibility purposes.

Some scripts assume:
- preprocessed signal matrices,
- local directory structures,
- and external genomics tools already installed.

Minor path modification may be required before execution.

---

# Citation

If you use this repository or PATTY in your work, please cite:

> Hu S. et al. PATTY: correction of open-chromatin bias in CUT&Tag profiling. (under review)

---

# Contact

Sheng’en Hu  
Department of Genome Sciences  
University of Virginia

GitHub Issues can also be used for questions or bug reports.
