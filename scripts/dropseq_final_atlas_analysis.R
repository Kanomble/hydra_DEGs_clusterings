library(Seurat)
# read the finalized atlas table

finalized_filepath <- "data/05_hydraAtlasReMap/04_finalize/nonDubLabeledSeurat.rds"
output_filepath <- "results/processedData/genes_to_cell_type_table.tsf"

atlas <- readRDS(finalized_filepath)
# counts <- atlas@assays$RNA@counts
norm_counts <- atlas@assays$SCT@data
ct <- atlas@active.ident
colnames(norm_counts) <- ct


# Get the unique column names
unique_col_names <- unique(colnames(norm_counts))

# Create initial dataframe
result <- data.frame()
col_name <- unique_col_names[1]
col_indices <- grep(col_name, colnames(norm_counts))
row_means <- rowMeans(norm_counts[, col_indices]) # rowMeans
df <- data.frame(row_means)
result <- merge(result, df, by = 0, all = TRUE)
rownames(result) <- result$Row.names
result$Row.names <- NULL
names(result)[names(result) == 'row_means'] <- col_name

# Testing potential cutoff values for assigning genes to celltypes 
# --> this is now done in a python script
#cutoff_values <- c(col_name)
#cutoff_values <- c(cutoff_values, max(unlist(result[col_name])))
#cutoff_values <- c(cutoff_values, mean(unlist(result[col_name])))
#cutoff_values <- c(cutoff_values, median(unlist(result[col_name])))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.001))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.005))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.01))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.02))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.05))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.1))
#cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.2))

# Loop over the unique column names and sum the values for each column
# Fill initial dataframe
for (col_name in unique_col_names[2:length(unique_col_names)]) {
  # Get the column indices with the current column name
  col_indices <- grep(col_name, colnames(norm_counts))
  
  # Sum the values of the columns with the current column name
  row_means <- rowMeans(norm_counts[, col_indices])
  df <- data.frame(row_means)
  
  # Add the row sums to the result data frame
  result <- merge(result, df, by = 0, all = TRUE)
  # Replace the column header "row_sums" by col_name
  names(result)[names(result) == 'row_means'] <- col_name
  # Reset the index column to contain gene names
  rownames(result) <- result$Row.names
  result$Row.names <- NULL
  
  # Different cutoff values will result in different dataframes ...
  #cutoff_values <- c(cutoff_values, col_name)
  #cutoff_values <- c(cutoff_values, max(unlist(result[col_name])))
  #cutoff_values <- c(cutoff_values, mean(unlist(result[col_name])))
  #cutoff_values <- c(cutoff_values, median(unlist(result[col_name])))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.001))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.005))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.01))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.02))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.05))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.1))
  #cutoff_values <- c(cutoff_values, sum(unlist(result[col_name]) <= 0.2))
  # hist(unlist(result[col_name]))
}
# Write result dataframe to disk
write.table(result, file=output_filepath, sep="\t", quote=F)