# A mixed-survival model of Lassa fever identifies correlates of outcome and severity

> Manuscript in preparation

## Methods

### Targeted transcriptomics (Nanostring)

Targeted transcriptomic analysis of the expression of 770 host mRNAs was quantified via the Nanostring NHP Immunology v2 panel (Bruker #531-115000276). The concentration and purity of whole blood RNA extracted from AVL following RT-qPCR were assessed using a Nanodrop spectrophotometer (Thermo Fisher Scientific). Samples were then diluted per the manufacturer’s instructions. Samples with <20 ng/µL RNA were used as-is. All samples were then prepared using the nCounter Prep Station and assayed via the nCounter Pro Analysis System (Bruker) according to the manufacturer’s instructions. Output files were loaded into [nSolver v4.0](https://nanostring.com/products/ncounter-analysis-system/ncounter-analysis-solutions/), and background thresholding was performed using the default parameters. Thresholded count matrices were exported from nSolver and analyzed with [limma v3.62.1](https://doi.org/10.1093/nar/gkv007) in R v4.4.2. Significantly upregulated or downregulated mRNAs were filtered by padj < 0.05 and log2 fold change > 1 or <-1.

### Bead-based multiplex assays
Circulating cytokines, thrombosis markers, and fibrinolysis analytes were measured by LEGENDplex™ bead-based immunoassays (BioLegend). Gamma-irradiated plasma samples were assessed in duplicate for each NHP Inflammation (#741491, 1:4 dilution), Human Thrombosis (#740892, 1:100 dilution), NHP chemokine/cytokine (#740388, 1:4 dilution), Human Immune Checkpoint (#740962, 1:2 dilution), Human Vascular Inflammation (#740590, 1:10,000) and Human Fibrinolysis (#740761, 1:40,000 dilution) panel according to the manufacturer instructions. Assay standards were mixed in batches and aliquoted across all plates to ensure consistency. Optional wash steps were incorporated to reduce the background signal. Assay samples were analyzed on an Accuri C6 Plus flow cytometer (BD Biosciences). Raw `.fcs` files from each assay were imported into LEGENDplex Qognit cloud-based Data Analysis Software suite (BioLegend), which automatically determined analyte concentrations in experimental samples following 5-parameter logistic regression curve fitting of assay standards. Analyte concentration data from assayed plasma samples were exported from Qognit and analyzed in R v4.3.1 using [limma v3.62.1](https://doi.org/10.1093/nar/gkv007). Significantly upregulated or downregulated proteins were filtered by padj < 0.05, and log2 fold change > 1 or <- 1.

## Data availability

The LV0043 genome consensus sequence was deposited in GenBank (`PX279630` and `PX279631`). Nanostring transcriptomic data were deposited in NCBI GEO (` TBD`). The unformatted output figures are in [`analysis`](analysis/). 
