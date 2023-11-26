---
editor_options: 
  markdown: 
    wrap: 72
---

## Project Contributors

Student ids and matching GitHub usernames for the group participants are
listed below:

-   s222845 - K0nuc1

-   s194562 - mathildewbech

-   s231488 - KKirdey

-   s215023 - lai1a

## Final project for course: R for Bio data Science

This project is built upon the scientific paper ***Bioinformatic
analysis reveals an exosomal miRNA-mRNA network in colorectal cancer***
written by Jun Ma et al. [1]

### Data retrieval - TCGA dataset

mRNA and miRNA data of colon adenocarcinoma and rectal adenocarcinoma
was gathered from the project TCGA-COAD from the GDC data portal:
<https://portal.gdc.cancer.gov/projects/TCGA-COAD> by using the R
library **TCGABiolinks**.

The two files ***analyte.tsv*** and ***clinical.tsv*** were first
retrieved from the "Biospecimen" archive on the specified webpage, and
further used to retrieve the read count datasets from the online data
portal. Both files should be downloaded and placed in the
**`"/TCGA_data/_raw/"`** data folder for the data loading to happen
successfully.

By running the script TCGA_load.qmd, all necessary data will be
downloaded and placed in the **`"/TCGA_data/_raw/"`** data folder and
three data and metadata files will be saved under `"/TCGA_data/"` for
further analysis.

### Data retrieval - GSE dataset

Exosomal miRNAs data were downloaded from the GEO dataset GSE39833 by
using the R library **GEOquery**. No files should be downloaded before
running the data loading step, as the library directly loads the data
into the R environment.

## Citations

[1] Ma, J., Wang, P., Huang, L. *et al.* Bioinformatic analysis reveals
an exosomal miRNA-mRNA network in colorectal cancer. *BMC Med
Genomics* **14**, 60 (2021).
<https://doi.org/10.1186/s12920-021-00905-2>
