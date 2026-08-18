# Spatial Transcriptomics Identifies a Nerve-Proximal Neuroimmune Program in HPV-Negative HNSCC

This repository contains code and supporting documentation for the manuscript:

**Spatial Transcriptomics Identifies a Nerve-Proximal Neuroimmune Program in HPV-Negative Head and Neck Squamous Cell Carcinoma**

**Authors:** Riya Chhabra, Alfred Kao, Reena Ding, Suravi Bajaj, Symphony Griffith Jackson, Wei Tse Li, Daniel J. John, Jessica Wang-Rodriguez, and Weg M. Ongkeko.

## Overview

Perineural invasion (PNI) is a clinically important feature of aggressive head and neck squamous cell carcinoma (HNSCC). This project integrates single-cell-resolution Xenium spatial transcriptomics, complementary Visium spatial transcriptomics, and TCGA-HNSC bulk transcriptomic survival analysis to characterize the tumor–nerve microenvironment, HPV-stratified nerve-proximal tumor behavior, immune organization, and the clinical relevance of a spatially derived NFE2L2/MDM2/PPARG gene program.

The survival analysis is intended as an exploratory clinical association rather than an externally validated prognostic biomarker.

## Analyses included

The analysis workflow includes:

1. Xenium preprocessing and quality control
2. Cell-type annotation and validation
3. Nerve-associated cell scoring and spatial zone assignment
4. Tumor–nerve proximity and exploratory perineural invasion index analyses
5. Differential expression and patient-adjusted pseudobulk analyses of nerve-proximal versus distal tumor cells
6. Visium spatial transcriptomic analysis, UCell PNI/EMT signature scoring, and SpatialCellChat ligand–receptor inference
7. TCGA-HNSC survival analysis of the spatially derived NFE2L2/MDM2/PPARG composite signature
8. Figure generation and supplementary analyses

## Public datasets

The study uses publicly available datasets:

- **Xenium HNSCC spatial transcriptomics:** GEO accession GSE300147
- **Visium HNSCC spatial transcriptomics:** GEO accessions GSE281978 and GSE181300
- **Bulk transcriptomic survival analysis:** TCGA-HNSC PanCancer Atlas 2018 via cBioPortal

Raw datasets should be downloaded from their original repositories. This repository provides analysis code and workflow documentation needed to reproduce the manuscript analyses.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── environment.yml
├── docs/
└── scripts/
```

## Reproducibility notes

Because the source datasets are large and hosted externally, raw data files are not included directly in this repository. Users should download the raw data from GEO and TCGA-HNSC/cBioPortal using the accessions listed above.

Before running analyses, users may need to create local directories for raw data, processed data, results, and figures, and adjust file paths in the scripts to match their local environment.

The Visium workflow was cleaned from the analysis scripts used to generate the reported results. Sample identities and SpatialCellChat inputs are made explicit in the repository scripts for reproducibility.

## Software environment

Python dependencies are listed in `requirements.txt` and `environment.yml`.

The Visium analyses additionally use R packages including Seurat, UCell, harmony, SpatialCellChat, Matrix, ggplot2, patchwork, and tidyverse. Exact package installation may depend on the user's R/Seurat environment.

## Citation

If you use this repository, please cite the associated manuscript and the archived Zenodo release.

Zenodo archive: **10.5281/zenodo.20953858**

## License

Code in this repository is released under the MIT License. Dataset access and reuse are governed by the terms of the original data repositories and source studies.
