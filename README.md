# Spatial Neuroimmune Crosstalk Driving Perineural Invasion in HNSCC

This repository contains code and supporting documentation for the manuscript:

**Spatial Neuroimmune Crosstalk Driving Perineural Invasion in Head and Neck Squamous Cell Carcinoma**

Authors: Riya Chhabra, Alfred Kao, Reena Ding, Suravi Bajaj, Symphony Griffith Jackson, Wei Tse Li, Daniel J. John, Jessica Wang-Rodriguez, and Weg M. Ongkeko.

## Overview

Perineural invasion (PNI) is a clinically important feature of aggressive head and neck squamous cell carcinoma (HNSCC). This project integrates single-cell-resolution Xenium spatial transcriptomics, Visium spatial transcriptomics, and bulk transcriptomic validation cohorts to characterize the tumor–nerve microenvironment, HPV-stratified nerve-proximal tumor behavior, immune organization, and prognostic relevance of a spatially derived NFE2L2/MDM2/PPARG gene signature.

## Analysis workflow

| Manuscript analysis | Script |
|---|---|
| Xenium preprocessing and annotation | scripts/01-10 |
| Visium preprocessing, UCell scoring, EMT analysis | scripts/11_visium_preprocessing_signature_scoring.R |
| SpatialCellChat analysis | scripts/12_spatialcellchat_visium.R |
| TCGA-HNSC survival analysis | scripts/13a_tcga_hnsc_survival_primarytumor.py |
| TCGA multivariable Cox analysis | scripts/13b_tcga_hnsc_multivariable_cox.py |
| GSE65858 external survival analysis | scripts/14_gse65858_survival.py |

## Public datasets

The study uses publicly available datasets:

* Xenium HNSCC spatial transcriptomics: GEO accession GSE300147
* Visium HNSCC spatial transcriptomics: GEO accessions GSE281978 and GSE181300
* Bulk transcriptomic validation: TCGA-HNSC PanCancer Atlas 2018 via cBioPortal
* Independent bulk transcriptomic validation: GEO accession GSE65858

Raw datasets should be downloaded from their original repositories. This repository is intended to provide analysis code and workflow documentation needed to reproduce the manuscript analyses.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .zenodo.json
├── requirements.txt
├── environment.yml
├── docs/
└── scripts/
```

## Reproducibility notes

Because the source datasets are large and hosted externally, raw data files are not included directly in this repository. Users should download the raw data from GEO, TCGA-HNSC/cBioPortal, and related public repositories using the accessions listed above.

Before running analyses, users may need to create local directories for raw data, processed data, results, and figures, and then adjust file paths in the scripts to match their local environment.

## Software environment

Python and R package requirements are listed in `requirements.txt` and `environment.yml`.

## Citation

If you use this repository, please cite the associated manuscript and the archived Zenodo release:

https://doi.org/10.5281/zenodo.20953858

## License

Code in this repository is released under the MIT License. Dataset access and reuse are governed by the terms of the original data repositories and source studies.
