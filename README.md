# PATTY Manuscript Source Data and Figure Generation Code

This repository contains the source data and analysis/plotting scripts used in the PATTY manuscript, mainly focusing on reproducing manuscript figures, supplementary figures, and associated source data files.

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

## `sourceData/`

Contains source data files corresponding to manuscript figures and supplementary figures.

These files were used to generate the figures shown in the manuscript.

---

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

# Notes

- Most scripts are figure-oriented and were used directly for manuscript figure generation.
- Data in many figures can be directly retrieved from sourceData/.
- Some large intermediate files and raw sequencing datasets are not included in this repository.

