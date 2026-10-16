# Hydra Cell-Type Clustering

This repository contains scripts for summarizing single-cell expression by cell type and assigning genes to cell-type clusters.

## Workflow

1. **Summarize expression by cell type.** The [dropseq_final_atlas_analysis.R script](https://github.com/Kanomble/hydra_DEGs_clusterings/blob/main/scripts/dropseq_final_atlas_analysis.R) reads the finalized Seurat atlas, uses normalized expression values, and calculates the mean expression of each gene across the cells in each cell type. It writes `genes_to_cell_type_table.tsf`.

2. **Assign genes to cell types.** The [create_sc_atlas_cluser_table.ipynb notebook](https://github.com/Kanomble/hydra_DEGs_clusterings/blob/main/scripts/create_sc_atlas_cluser_table.ipynb) uses the cell-type mean-expression table to identify genes whose expression in a cell type exceeds their mean across cell types, and formats the assignments as a gene-to-cluster table.

3. **Build clusters and identify unique assignments.** The [get_cluster.R script](https://github.com/Kanomble/hydra_DEGs_clusterings/blob/main/scripts/get_cluster.R) combines cell types into clusters and higher-level categories, then creates clustering tables that distinguish unique from non-unique gene assignments. Its outputs include `hvaep_clustering_table.tsv` and `hvaep_clustering_table_unique.tsv`.

4. **Filter for unique cell-type assignments.** The `get_unique_cell_types.R` script applies the selected DESeq2 comparisons and filters to the clustered gene annotations. Within each hierarchy block, genes assigned to more than two categories are excluded; the resulting table is written to `results/processedData/unique_celltypes.table`.

## Key outputs

- `genes_to_cell_type_table.tsf`: mean normalized expression by gene and cell type.
- `hvaep_cell_type_to_gene_cluster_table.tsf`: gene-to-cell-type assignments.
- `hvaep_clustering_table.tsv` and `hvaep_clustering_table_unique.tsv`: cluster assignments and unique assignments.
- `unique_celltypes.table`: cell-type annotations filtered using the selected DESeq2 comparisons and hierarchy-specific assignment threshold.
