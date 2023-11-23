library(edgeR)
library(statmod)

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