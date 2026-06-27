# HNSCC Perineural Invasion Spatial Transcriptomics

This repository contains analysis code for the manuscript:

**Spatial Neuroimmune Crosstalk Driving Perineural Invasion in Head and Neck Squamous Cell Carcinoma**

The study integrates single-cell-resolution Xenium spatial transcriptomics, Visium spatial transcriptomics, and bulk transcriptomic validation cohorts to characterize the tumor–nerve microenvironment in head and neck squamous cell carcinoma (HNSCC). Analyses include nerve-associated cell scoring, spatial zone assignment, tumor–nerve proximity analysis, perineural invasion index calculation, differential expression analysis, EMT and PNI signature scoring, inferred ligand–receptor analysis, and survival validation of a spatially derived NFE2L2/MDM2/PPARG signature.

## Data availability

No raw sequencing or spatial transcriptomic data are stored in this repository. All datasets analyzed in the manuscript are publicly available from the following sources:

* Xenium spatial transcriptomic data: GEO accession GSE300147
* Visium spatial transcriptomic data: GEO accessions GSE281978 and GSE181300
* Bulk transcriptomic validation data: TCGA-HNSC via cBioPortal, PanCancer Atlas 2018
* Independent bulk validation cohort: GEO accession GSE65858

## Repository structure

```text
hnscc-pni-spatial-transcriptomics/
├── README.md
├── LICENSE
├── CITATION.cff
├── .zenodo.json
├── requirements.txt
├── environment.yml
└── scripts/
    ├── 01_xenium_qc_preprocessing_GSE300147.py
    ├── 02_nerve_identification_PNI_zones.py
    ├── 03_gene_signatures_DEG.py
    └── additional analysis scripts
```

## Code overview

The scripts in this repository support the following analyses:

1. Xenium quality control, preprocessing, normalization, integration, clustering, and cell-type annotation.
2. Identification of Schwann/neural marker-expressing nerve-associated cells.
3. Spatial assignment of cells into nerve-associated, perineural, peritumoral/local microenvironmental, and distal zones.
4. Quantification of tumor–nerve proximity and calculation of a per-patient exploratory PNI index.
5. Differential expression analysis of perineural versus distal tumor cells.
6. Visium-based EMT and PNI signature scoring.
7. SpatialCellChat ligand–receptor inference.
8. TCGA-HNSC and GSE65858 survival validation of the NFE2L2/MDM2/PPARG signature.

## Software requirements

Python and R package requirements are listed in `requirements.txt` and `environment.yml`. Package versions reflect the computational environment used for the analyses where available.

## Citation

If you use this code, please cite the associated manuscript and archived Zenodo release:

https://doi.org/10.5281/zenodo.20953858

## Contact

For questions about the manuscript or repository, please contact the corresponding author listed in the manuscript.
