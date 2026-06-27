# Spatial Neuroimmune Crosstalk Driving Perineural Invasion in HNSCC

This repository contains code and supporting materials for the manuscript:

**Spatial Neuroimmune Crosstalk Driving Perineural Invasion in Head and Neck Squamous Cell Carcinoma**

Authors: Riya Chhabra, Alfred Kao, Reena Ding, Suravi Bajaj, Symphony Griffith Jackson, Wei Tse Li, Daniel J. John, Jessica Wang-Rodriguez, and Weg M. Ongkeko.

## Overview

Perineural invasion (PNI) is a clinically important feature of aggressive head and neck squamous cell carcinoma (HNSCC). This project integrates single-cell-resolution Xenium spatial transcriptomics, Visium spatial transcriptomics, and bulk transcriptomic validation cohorts to characterize the tumor–nerve microenvironment, HPV-stratified nerve-proximal tumor behavior, immune organization, and prognostic relevance of a spatially derived gene signature.

## Analyses included

The analysis workflow is organized around the following components:

1. Xenium preprocessing and quality control
2. Cell type annotation and validation
3. Nerve-associated cell scoring and spatial zone assignment
4. Tumor–nerve proximity and perineural invasion index analyses
5. Differential expression and pseudobulk analyses of perineural versus distal tumor cells
6. Visium preprocessing, UCell signature scoring, EMT/PNI enrichment analyses, and SpatialCellChat inference
7. TCGA-HNSC and GSE65858 bulk transcriptomic survival validation of the NFE2L2/MDM2/PPARG signature
8. Figure generation

## Public datasets

The study uses publicly available datasets:

- Xenium HNSCC spatial transcriptomics: GEO accession **GSE300147**
- Visium HNSCC spatial transcriptomics: GEO accessions **GSE281978** and **GSE181300**
- Bulk transcriptomic validation: **TCGA-HNSC PanCancer Atlas 2018** via cBioPortal
- Independent bulk transcriptomic validation: GEO accession **GSE65858**

Raw datasets should be downloaded from their original repositories. This repository is intended to provide analysis code, workflow documentation, and derived non-identifiable outputs needed to reproduce the manuscript analyses.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .zenodo.json
├── .gitignore
├── requirements.txt
├── environment.yml
├── data/
│   └── README.md
├── scripts/
│   └── README.md
├── figures/
│   └── README.md
├── results/
│   └── README.md
└── docs/
    └── manuscript_data_availability_statement.md
```

## Reproducibility notes

Because several source datasets are large and hosted externally, raw data files are not included directly in this repository. Users should download the raw data from GEO, TCGA/cBioPortal, and related public repositories using the accessions listed above.

Before running analyses, create the expected local directories:

```bash
mkdir -p data/raw data/processed results figures
```

Then place downloaded raw files under `data/raw/` and run the scripts in the order described in `scripts/README.md`.

## Software environment

Example Python and R package requirements are listed in `requirements.txt` and `environment.yml`. Exact package versions should be updated to match the final computational environment used for the manuscript.

## Citation

If you use this repository, please cite the associated manuscript and the archived Zenodo DOI once available.

## License

Code in this repository is released under the MIT License. Dataset access and reuse are governed by the terms of the original data repositories and source studies.
