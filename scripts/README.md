# Scripts

This folder contains the analysis scripts used for the manuscript:

**Spatial Neuroimmune Crosstalk Driving Perineural Invasion in Head and Neck Squamous Cell Carcinoma**

The scripts support Xenium spatial transcriptomic preprocessing, nerve-associated cell identification, spatial zone assignment, differential expression analysis, cell-type annotation and validation, pseudobulk analysis, neighborhood analysis, dendritic cell validation, and figure generation.

## Script overview

### `01_xenium_qc_preprocessing_GSE300147.py`

Performs quality control, preprocessing, normalization, dimensionality reduction, clustering, and preparation of the Xenium spatial transcriptomic dataset from GSE300147.

### `02_nerve_identification_PNI_zones.py`

Identifies Schwann/neural marker-expressing nerve-associated cells, calculates distance to the nearest nerve-associated cell, and assigns cells to nerve-associated, perineural, peritumoral/local microenvironmental, or distal spatial zones.

### `02b_nerve_threshold_justification.py`

Evaluates nerve-associated marker score thresholds and supports the threshold sensitivity analyses used to justify nerve-associated cell classification.

### `03_gene_signatures_DEG.py`

Performs gene signature scoring and differential expression analysis of perineural versus distal tumor cells, including analyses related to the spatially enriched perineural tumor-cell program.

### `04_cell_type_annotation.py`

Annotates major tumor, immune, stromal, endothelial, and mast cell populations using canonical marker genes and supports downstream cell-type-specific analyses.

### `05_additional_figures.py`

Generates additional manuscript and supplementary figures related to spatial organization, cell-type composition, and downstream analyses.

### `05b_fix_figures.py`

Applies final figure formatting and plotting adjustments for manuscript-ready visualizations.

### `06_puram_validation_v3.py`

Performs reference-based validation of cell-type annotations using the published HNSCC single-cell reference dataset from Puram et al.

### `07_replicate_concordance.py`

Assesses concordance across samples, replicates, or related analysis outputs to support robustness of observed spatial and transcriptional patterns.

### `08_pseudobulk_final.py`

Performs patient-level pseudobulk differential expression analysis to account for inter-patient variability in perineural versus distal tumor-cell comparisons.

### `08b_pseudobulk_HPVneg.py`

Performs HPV-negative-stratified pseudobulk analysis to evaluate nerve-proximal transcriptional programs in HPV-negative disease.

### `09_neighbourhood_analysis.py`

Performs spatial neighborhood enrichment analysis around nerve-associated cells and evaluates cell-type enrichment or depletion in the nerve-proximal microenvironment.

### `10_lamp3_dc_validation.py`

Validates the DC_LAMP3/mature dendritic cell annotation using curated marker gene scoring and marker specificity analysis.

### 11_visium_preprocessing_signature_scoring.R
Processes seven public HNSCC Visium tissue sections using Seurat,
including QC, SCTransform normalization, Harmony integration,
cell-type annotation, UCell-based nerve/PNI scoring, EMT scoring,
and generation of Figure 5A-D analyses.

### 12_spatialcellchat_visium.R
Performs SpatialCellChat analysis on three Visium sections with
sufficient nerve-associated spots (GSM8633893, GSM5494475,
GSM5494476), including spatially constrained ligand-receptor
inference and Figure 5E visualization.

### 13a_tcga_hnsc_survival_primarytumor.py
Reproduces the TCGA-HNSC primary-tumor survival analysis for the
NFE2L2/MDM2/PPARG composite signature, including HPV-stratified
Kaplan-Meier analyses and REMARK cohort/event reporting.

### 13b_tcga_hnsc_multivariable_cox.py
Performs the primary multivariable TCGA-HNSC Cox regression
adjusting for age and AJCC stage using complete-case analysis.

### 14_gse65858_survival.py
Reproduces the GSE65858 external survival analysis of the
NFE2L2/MDM2/PPARG signature, including HPV-stratified 5-year
Kaplan-Meier analyses and REMARK cohort/event reporting.

## Notes

Raw spatial transcriptomic and bulk transcriptomic datasets are not included in this repository. Users should download the required public datasets from GEO, TCGA-HNSC/cBioPortal, and GSE65858 as described in the main repository README and manuscript.

The scripts are intended to document and reproduce the analyses reported in the manuscript. Some file paths may need to be adjusted depending on the user’s local directory structure.
