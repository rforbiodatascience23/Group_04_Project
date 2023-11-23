library(edgeR)
library(statmod)

# Functions used in TCGA_load
# Function to the IDs of specific type (cancer and normal)

get_sample_ids <- function(RNA_samples_data, sample_type, num_samples) {
  
  filtered_samples <- RNA_samples_data |>
    filter(type == sample_type) |>
    head(num_samples)
  
  id_cancer_patients <- filtered_samples$analyte_submitter_id
  
  return(id_cancer_patients)
}

# Function to retrieve and prepare data from TCGA

retrieve_and_prepare_data <- function(
    project,
    data_category,
    data_type,
    workflow_type,
    id_cancer_patients,
    directory_prefix,
    use_prepare = TRUE
) {
  # Query to specify the data to get
  query <- GDCquery(
    project = project,
    data.category = data_category,
    data.type = data_type,
    workflow.type = workflow_type,
    barcode = id_cancer_patients
  )
  
  # Downloading the samples
  GDCdownload(
    query = query,
    method = "api",
    directory = paste0(directory_prefix, "_", project),
    files.per.chunk = 50
  )
  
  # Preparing data in case of miRNAs
  if (use_prepare) {
    data <- GDCprepare(query, directory = paste0(directory_prefix, "_", project))
  } else {
    data <- NULL
  }
  return(data)
}


# Functions used in TCGA_augment 
# Function to create normalized expression datasets.

process_data <- function(data) {
  # Extract IDs and expression data
  ids <- data[, 1]  # IDs are in the first column
  expression_data <- data[, -1]
  
  # Imputation to estimate missing values based on the mean of available data
  expression_data_imputed <- apply(expression_data, 2, 
                                   function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  )
  
  # Create an edgeR object
  y <- DGEList(counts = expression_data_imputed, genes = ids)
  
  # Normalization step
  y <- calcNormFactors(y)
  
  # Return the result
  return(y)
}


# Function to run edgeR differential expression analysis.

run_edgeR_analysis <- function(data_y,
                               design_formula,
                               metadata, 
                               feature_data) {
  # Create model matrix
  # design formula is the same in this case, but in case of other names, it can be changed
  design <- model.matrix(as.formula(design_formula), 
                         data = metadata)
  
  # Set row names
  rownames(design) <- data_y$samples[, 0]
  
  # Estimate dispersion
  data_y <- estimateDisp(data_y, design)
  
  # Fit model
  fit <- glmQLFit(data_y, design, robust = TRUE)
  
  # Likelihood ratio test
  qlt <- glmQLFTest(fit)
  
  # Getting top genes and saving table in the variable. 
  # The dimension of the original dataset is applied.
  topgenes <- topTags(qlt, n = dim(feature_data)[[1]])
  
  return(topgenes)
}

# Function to pivot the datasets for the further merging with metadata.

pivot_and_filter <- function(data, 
                             id_column, 
                             top_genes_list) {
  filtered_data <- data |>
    filter({{ id_column }} %in% top_genes_list) |>
    pivot_longer(
      cols = -{{ id_column }},
      names_to = "TCGA_ID",
      values_to = "log_reads"
    )
  
  return(filtered_data)
}