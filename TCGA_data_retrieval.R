setwd("/Users/kseniakirdey/Desktop/DTU/R_for_bio_data_science/lab_10_Project")

library(TCGAbiolinks)
library(tidyverse)
library(stringr)

analyte <- read.delim("analyte.tsv")
#listSamples <- head(analyte$sample_submitter_id, 10)

RNA_samples <- subset(analyte, analyte_type == "RNA")

RNA_samples <- RNA_samples |>
  mutate(tissue_type = str_sub(sample_submitter_id, -3, -3)) |>
  mutate(
    type = case_when(tissue_type == 0 ~ "cancer",
                     (tissue_type == 1 | tissue_type == 2) ~ "normal",
                     .default = "other")
  ) |>
  arrange(type)
  
cancer_RNA_samples <- RNA_samples |> 
  filter(type == "cancer") |>
  head(10)

normal_RNA_samples <- RNA_samples |> 
  filter(type == "normal") |>
  head(10)

id_cancer_patients_normal <- normal_RNA_samples$analyte_submitter_id
id_cancer_patients_cancer <- cancer_RNA_samples$analyte_submitter_id

# MiRNA data

query_miRNA_cancer <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  barcode = id_cancer_patients_normal
)

GDCdownload(
  query = query_miRNA_cancer,
  method = "api",
  directory = "samples_miRNA_COAD_cancer",
  files.per.chunk = 50
)

# Prepare miRNA expression matrix
miRNA_data <- GDCprepare(query_miRNA, directory = "samples_miRNA_COAD")

# MRNA data
# Query for mRNA data using the same samples
query_mRNA <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = id_cancer_patients
)

# Download mRNA data (if not already downloaded)
GDCdownload(
  query = query_mRNA,
  method = "api",
  directory = "samples_mRNA_COAD",
  files.per.chunk = 50
)

mRNA_data <- GDCprepare(query_mRNA, directory = "samples_mRNA_COAD")