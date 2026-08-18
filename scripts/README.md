# Scripts

This folder contains the analysis scripts used for the manuscript:

**Spatial Transcriptomics Identifies a Nerve-Proximal Neuroimmune Program in HPV-Negative Head and Neck Squamous Cell Carcinoma**

The scripts support Xenium spatial transcriptomic preprocessing, nerve-associated cell identification, spatial zone assignment, differential expression analysis, cell-type annotation and validation, patient-adjusted pseudobulk analysis, spatial neighborhood analysis, dendritic-cell validation, complementary Visium analysis, SpatialCellChat ligand–receptor inference, TCGA-HNSC survival analysis, and figure generation.

## Script overview

### `01_xenium_qc_preprocessing_GSE300147.py`

Performs quality control, preprocessing, normalization, dimensionality reduction, clustering, and preparation of the Xenium spatial transcriptomic dataset from GSE300147.

### `02_nerve_identification_PNI_zones.py`

Identifies Schwann/neural marker-expressing nerve-associated cells, calculates distance to the nearest nerve-associated cell, and assigns cells to nerve-associated, perineural, peritumoral/local microenvironmental, or distal spatial zones.

### `02b_nerve_threshold_justification.py`

Evaluates nerve-associated marker-score thresholds and supports the threshold-sensitivity analyses used to justify nerve-associated cell classification.

### `03_gene_signatures_DEG.py`

Performs gene-signature scoring and differential expression analysis of nerve-proximal versus distal tumor cells, including analyses of the spatially enriched nerve-proximal tumor-cell program.

### `04_cell_type_annotation.py`

Annotates major tumor, immune, stromal, endothelial, and mast-cell populations using canonical marker genes and supports downstream cell-type-specific analyses.

### `05_additional_figures.py`

Generates additional manuscript and supplementary figures related to spatial organization, cell-type composition, and downstream analyses.

### `05b_fix_figures.py`

Applies final figure-formatting and plotting adjustments for manuscript-ready visualizations.

### `06_puram_validation_v3.py`

Performs reference-based validation of cell-type annotations using the published HNSCC single-cell reference dataset from Puram et al.

### `07_replicate_concordance.py`

Assesses concordance across samples, replicates, or related analysis outputs to support robustness of observed spatial and transcriptional patterns.

### `08_pseudobulk_final.py`

Performs patient-level pseudobulk differential expression analysis to account for inter-patient variability in nerve-proximal versus distal tumor-cell comparisons.

### `08b_pseudobulk_HPVneg.py`

Performs HPV-negative-stratified pseudobulk analysis to evaluate nerve-proximal transcriptional programs in HPV-negative disease.

### `09_neighbourhood_analysis.py`

Performs spatial neighborhood enrichment analysis around nerve-associated cells and evaluates cell-type enrichment or depletion in the nerve-proximal microenvironment.

### `10_lamp3_dc_validation.py`

Validates the DC_LAMP3/mature-migratory dendritic-cell annotation using curated marker-gene scoring and marker-specificity analysis.

### `11_visium_preprocessing_ucell.R`

Processes seven public HNSCC Visium tissue sections using Seurat, including quality control, SCTransform normalization, Harmony integration, cell-type annotation, UCell-based nerve/PNI scoring, EMT scoring, and analyses supporting Figure 5A–D.

### `12_visium_spatialcellchat.R`

Performs SpatialCellChat analysis on three Visium sections with sufficient nerve-associated spots (GSM8633893, GSM5494475, and GSM5494476), including spatially constrained ligand–receptor inference and analyses supporting Figure 5E.

### `13a_tcga_hnsc_survival_primarytumor.py`

Reproduces the corrected TCGA-HNSC primary-tumor survival analysis for the NFE2L2/MDM2/PPARG composite signature. RNA-seq data are restricted to primary solid tumor samples before patient-level barcode collapsing, followed by HPV-stratified Kaplan–Meier analysis and cohort/event reporting.

### `13b_tcga_hnsc_multivariable_cox.py`

Performs the primary multivariable TCGA-HNSC Cox regression analysis of the continuous NFE2L2/MDM2/PPARG composite score in HPV-negative disease, adjusting for age and categorical AJCC stage using complete-case analysis.

## Notes

Raw spatial transcriptomic and bulk transcriptomic datasets are not included in this repository. Users should download the required public datasets from GEO and TCGA-HNSC/cBioPortal as described in the main repository README and manuscript.

The scripts are intended to document and reproduce the analyses reported in the manuscript. Some file paths may need to be adjusted depending on the user's local directory structure. The Visium scripts document the analysis workflow used for the reported results; the underlying Visium data remain available from the original GEO repositories rather than being redistributed here.
