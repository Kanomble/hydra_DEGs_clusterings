# ============================================================
# DESeq2 + combined transcriptional-response PCA / PCoA
#
# Tailored to:
#   data/rsem_counts/sample_descriptions_pcoa.csv
#   data/rsem_counts/sample_controls_pcoa.csv
#
# Combined analysis:
#   - each point = one DESeq2 treatment-vs-control contrast
#   - each variable = one gene's log2 fold change
#   - PCA on the log2FC matrix
#   - Pearson-correlation heatmap
#   - PCoA using distance = 1 - Pearson correlation
#
# IMPORTANT:
# HydraRecolonization batch is confounded with treatment after
# the specified bad samples are removed, so the DESeq2 design is
# intentionally ~ condition, not ~ batch + condition.
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "tximport",
  "DESeq2",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
})


# ============================================================
# Input / output
# ============================================================

sample_file <- "data/rsem_counts/sample_descriptions_pcoa.csv"
controls_file <- "data/rsem_counts/sample_controls_pcoa.csv"

input_dir <- "data/rsem_counts"
output_dir <- "results/deseq2_rsem"

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dds_dir <- file.path(output_dir, "dds_objects")
combined_dir <- file.path(output_dir, "combined_analysis")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Settings
# ============================================================

# Keep genes with at least this many total counts within an experiment.
count_threshold <- 20

# For the combined contrast-level PCA/PCoA:
#   0   = use all genes with finite log2FC in every contrast
#   0.5 = keep genes with |log2FC| >= 0.5 in at least one contrast
combined_min_abs_log2fc <- 0

# Samples previously identified as bad in HydraRecolonization.
bad_recolonization_samples <- c(
  "I27673-S1",
  "I27674-S1",
  "I27675-S1",
  "I27678-S1",
  "I27681-S1",
  "I27672-S1"
)


# ============================================================
# Helper functions
# ============================================================

trim_character_columns <- function(df) {
  char_cols <- vapply(df, is.character, logical(1))
  df[char_cols] <- lapply(df[char_cols], trimws)
  df
}


safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}


save_sample_pca <- function(dds, experiment_name, output_file) {

  vsd <- varianceStabilizingTransformation(
    dds,
    blind = TRUE
  )

  pca_data <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
  )

  percent_var <- round(
    100 * attr(pca_data, "percentVar"),
    1
  )

  p <- ggplot(
    pca_data,
    aes(
      x = PC1,
      y = PC2,
      color = condition,
      label = name
    )
  ) +
    geom_point(size = 3) +
    geom_text(
      vjust = -0.8,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    labs(
      title = paste0(experiment_name, ": sample-level PCA"),
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Condition"
    ) +
    theme_classic()

  ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = 5.5,
    dpi = 150
  )

  invisible(p)
}


save_volcano_plot <- function(
    res_df,
    experiment_name,
    treatment_group,
    control_group,
    output_file
) {

  plot_df <- res_df[
    is.finite(res_df$log2FoldChange) &
      !is.na(res_df$padj) &
      is.finite(res_df$padj),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0) {
    warning(
      "No finite log2FoldChange/padj values for volcano plot: ",
      experiment_name, " ",
      treatment_group, " vs ", control_group
    )
    return(invisible(NULL))
  }

  plot_df$minus_log10_padj <- -log10(
    pmax(plot_df$padj, .Machine$double.xmin)
  )

  plot_df$significance <- "Not significant"
  plot_df$significance[
    plot_df$padj < 0.05 &
      abs(plot_df$log2FoldChange) >= 1
  ] <- "padj < 0.05 & |log2FC| >= 1"

  label_df <- plot_df[
    order(plot_df$padj, na.last = TRUE),
    ,
    drop = FALSE
  ]

  label_df <- head(label_df, 10)

  p <- ggplot(
    plot_df,
    aes(
      x = log2FoldChange,
      y = minus_log10_padj,
      color = significance
    )
  ) +
    geom_point(alpha = 0.65, size = 1.5) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_text(
      data = label_df,
      aes(label = gene),
      vjust = -0.6,
      size = 2.8,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    labs(
      title = paste0(
        experiment_name,
        ": ",
        treatment_group,
        " vs ",
        control_group
      ),
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_classic()

  ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = 5.5,
    dpi = 150
  )

  invisible(p)
}


# ============================================================
# Read metadata
# ============================================================

samples <- read.csv(
  sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

controls <- read.csv(
  controls_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

samples <- trim_character_columns(samples)
controls <- trim_character_columns(controls)


# ============================================================
# Normalize control-file column names
# ============================================================

if ("treatment group" %in% names(controls)) {
  names(controls)[names(controls) == "treatment group"] <- "treatment_group"
}

if ("control group" %in% names(controls)) {
  names(controls)[names(controls) == "control group"] <- "control_group"
}


# ============================================================
# Validate metadata columns
# ============================================================

required_sample_columns <- c(
  "EXP",
  "ID",
  "sample",
  "treatment",
  "batch"
)

required_control_columns <- c(
  "EXP",
  "treatment_group",
  "control_group"
)

missing_sample_columns <- setdiff(
  required_sample_columns,
  names(samples)
)

missing_control_columns <- setdiff(
  required_control_columns,
  names(controls)
)

if (length(missing_sample_columns) > 0) {
  stop(
    "Missing columns in sample description file: ",
    paste(missing_sample_columns, collapse = ", ")
  )
}

if (length(missing_control_columns) > 0) {
  stop(
    "Missing columns in controls file: ",
    paste(missing_control_columns, collapse = ", ")
  )
}


# ============================================================
# Validate contrast definitions against sample metadata
# ============================================================

for (i in seq_len(nrow(controls))) {

  exp_name <- controls$EXP[i]
  treatment_group <- controls$treatment_group[i]
  control_group <- controls$control_group[i]

  if (!exp_name %in% samples$EXP) {
    stop(
      "Experiment in controls file is absent from sample file: ",
      exp_name
    )
  }

  exp_treatments <- unique(
    samples$treatment[samples$EXP == exp_name]
  )

  if (!treatment_group %in% exp_treatments) {
    stop(
      "Treatment '", treatment_group,
      "' is not present in experiment '", exp_name, "'."
    )
  }

  if (!control_group %in% exp_treatments) {
    stop(
      "Control '", control_group,
      "' is not present in experiment '", exp_name, "'."
    )
  }
}


# ============================================================
# Containers for results
# ============================================================

experiments <- unique(samples$EXP)

dds_objects <- list()
contrast_results <- list()
contrast_metadata <- list()


# ============================================================
# Main DESeq2 loop
# ============================================================

for (exp_name in experiments) {

  cat("\n")
  cat("============================================================\n")
  cat("EXPERIMENT: ", exp_name, "\n", sep = "")
  cat("============================================================\n")

  exp_samples <- samples[
    samples$EXP == exp_name,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Remove specified bad HydraRecolonization samples
  # ----------------------------------------------------------

  if (exp_name == "HydraRecolonization") {

    bad_present <- intersect(
      bad_recolonization_samples,
      exp_samples$ID
    )

    if (length(bad_present) > 0) {
      cat(
        "[*] Removing bad samples: ",
        paste(bad_present, collapse = ", "),
        "\n",
        sep = ""
      )

      exp_samples <- exp_samples[
        !exp_samples$ID %in% bad_present,
        ,
        drop = FALSE
      ]
    }
  }


  # ----------------------------------------------------------
  # Print treatment x batch structure
  # ----------------------------------------------------------

  cat("[*] Treatment x batch table:\n")
  print(
    with(
      exp_samples,
      table(treatment, batch)
    )
  )


  # ----------------------------------------------------------
  # RSEM input files
  # ----------------------------------------------------------

  files <- file.path(
    input_dir,
    exp_name,
    paste0(exp_samples$ID, ".genes.results")
  )

  names(files) <- exp_samples$ID

  missing_files <- files[
    !file.exists(files)
  ]

  if (length(missing_files) > 0) {
    stop(
      "Missing RSEM files for ", exp_name, ":\n",
      paste(missing_files, collapse = "\n")
    )
  }

  cat(
    "[*] Importing ",
    length(files),
    " RSEM files\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # tximport
  # ----------------------------------------------------------

  txi_rsem <- tximport(
    files,
    type = "rsem",
    txIn = FALSE,
    txOut = FALSE
  )


  # ----------------------------------------------------------
  # Sample table
  # ----------------------------------------------------------

  sample_table <- data.frame(
    condition = factor(exp_samples$treatment),
    batch = factor(exp_samples$batch),
    row.names = exp_samples$ID
  )

  sample_table <- sample_table[
    colnames(txi_rsem$counts),
    ,
    drop = FALSE
  ]

  if (!identical(
    rownames(sample_table),
    colnames(txi_rsem$counts)
  )) {
    stop(
      "Sample metadata and tximport columns do not match for ",
      exp_name
    )
  }


  # ----------------------------------------------------------
  # DESeq2 dataset
  #
  # Design is intentionally ~ condition.
  #
  # HydraAHL has only one batch.
  # HydraRecolonization has treatment/batch confounding after
  # the bad samples are removed, so batch cannot be estimated
  # independently from treatment.
  # ----------------------------------------------------------

  dds <- DESeqDataSetFromTximport(
    txi_rsem,
    colData = sample_table,
    design = ~ condition
  )


  # ----------------------------------------------------------
  # Count filter
  # ----------------------------------------------------------

  cat(
    "[*] Genes before count filter: ",
    nrow(dds),
    "\n",
    sep = ""
  )

  keep <- rowSums(counts(dds)) >= count_threshold

  dds <- dds[
    keep,
    ,
    drop = FALSE
  ]

  cat(
    "[*] Genes after count filter: ",
    nrow(dds),
    "\n",
    sep = ""
  )

  if (nrow(dds) == 0) {
    stop(
      "No genes remain after filtering in experiment ",
      exp_name
    )
  }

  dds$condition <- droplevels(dds$condition)
  dds$batch <- droplevels(dds$batch)


  # ----------------------------------------------------------
  # Run DESeq2
  # ----------------------------------------------------------

  cat("[*] Running DESeq2\n")

  dds <- DESeq(dds)

  dds_objects[[exp_name]] <- dds

  saveRDS(
    dds,
    file = file.path(
      dds_dir,
      paste0(
        safe_filename(exp_name),
        "_dds.rds"
      )
    )
  )


  # ----------------------------------------------------------
  # Sample-level PCA
  # ----------------------------------------------------------

  pca_file <- file.path(
    figures_dir,
    paste0(
      safe_filename(exp_name),
      "_sample_PCA.png"
    )
  )

  cat(
    "[*] Writing sample PCA: ",
    pca_file,
    "\n",
    sep = ""
  )

  save_sample_pca(
    dds = dds,
    experiment_name = exp_name,
    output_file = pca_file
  )


  # ==========================================================
  # Contrasts for this experiment
  # ==========================================================

  exp_controls <- controls[
    controls$EXP == exp_name,
    ,
    drop = FALSE
  ]

  if (nrow(exp_controls) == 0) {
    warning(
      "No contrasts defined for experiment ",
      exp_name
    )
    next
  }

  cat(
    "[*] Number of contrasts: ",
    nrow(exp_controls),
    "\n",
    sep = ""
  )


  for (contrast_i in seq_len(nrow(exp_controls))) {

    treatment_group <- exp_controls$treatment_group[contrast_i]
    control_group <- exp_controls$control_group[contrast_i]

    cat(
      "\n[*] Contrast: ",
      treatment_group,
      " vs ",
      control_group,
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # Verify contrast levels after sample removal
    # --------------------------------------------------------

    condition_levels <- levels(dds$condition)

    if (!treatment_group %in% condition_levels) {
      stop(
        "Treatment '", treatment_group,
        "' is absent after sample filtering in ",
        exp_name
      )
    }

    if (!control_group %in% condition_levels) {
      stop(
        "Control '", control_group,
        "' is absent after sample filtering in ",
        exp_name
      )
    }


    # --------------------------------------------------------
    # DESeq2 result
    # --------------------------------------------------------

    res <- results(
      dds,
      contrast = c(
        "condition",
        treatment_group,
        control_group
      )
    )

    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    # Put gene ID first in output table.
    res_df <- res_df[
      ,
      c(
        "gene",
        setdiff(names(res_df), "gene")
      ),
      drop = FALSE
    ]

    res_ordered <- res_df[
      order(res_df$pvalue, na.last = TRUE),
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # Contrast name
    # --------------------------------------------------------

    contrast_name <- paste0(
      exp_name,
      "__",
      treatment_group,
      "_vs_",
      control_group
    )

    if (contrast_name %in% names(contrast_results)) {
      stop(
        "Duplicate contrast name detected: ",
        contrast_name
      )
    }


    # --------------------------------------------------------
    # Store for combined analysis
    # --------------------------------------------------------

    contrast_results[[contrast_name]] <- res_df

    contrast_metadata[[contrast_name]] <- data.frame(
      contrast = contrast_name,
      experiment = exp_name,
      treatment = treatment_group,
      control = control_group,
      stringsAsFactors = FALSE
    )


    # --------------------------------------------------------
    # Save DESeq2 table
    # --------------------------------------------------------

    result_file <- file.path(
      tables_dir,
      paste0(
        safe_filename(exp_name),
        "_",
        safe_filename(treatment_group),
        "_vs_",
        safe_filename(control_group),
        "_results.csv"
      )
    )

    write.csv(
      res_ordered,
      file = result_file,
      row.names = FALSE
    )

    cat(
      "[*] Wrote DESeq2 table: ",
      result_file,
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # Volcano plot
    # --------------------------------------------------------

    volcano_file <- file.path(
      figures_dir,
      paste0(
        safe_filename(exp_name),
        "_",
        safe_filename(treatment_group),
        "_vs_",
        safe_filename(control_group),
        "_volcano.png"
      )
    )

    save_volcano_plot(
      res_df = res_df,
      experiment_name = exp_name,
      treatment_group = treatment_group,
      control_group = control_group,
      output_file = volcano_file
    )

    cat(
      "[*] Wrote volcano plot: ",
      volcano_file,
      "\n",
      sep = ""
    )
  }

  cat(
    "\n[*] Finished experiment: ",
    exp_name,
    "\n",
    sep = ""
  )
}

# ============================================================
# COMBINED SAMPLE-LEVEL PCoA
#
# One PCoA containing individual samples from:
#   - HydraAHL
#   - HydraRecolonization
#
# Each point = one individual RNA-seq sample
#
# Distance = 1 - Pearson correlation of VST expression
# ============================================================

cat("\n")
cat("============================================================\n")
cat("COMBINED SAMPLE-LEVEL PCoA\n")
cat("HydraAHL + HydraRecolonization\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Start from complete sample metadata
# ------------------------------------------------------------

all_samples <- samples


# ------------------------------------------------------------
# Remove the bad HydraRecolonization samples
# ------------------------------------------------------------

all_samples <- all_samples[
  !all_samples$ID %in% bad_recolonization_samples,
  ,
  drop = FALSE
]


cat(
  "[*] Number of individual samples: ",
  nrow(all_samples),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Show samples included
# ------------------------------------------------------------

print(
  with(
    all_samples,
    table(EXP, treatment)
  )
)


# ------------------------------------------------------------
# Generate RSEM file paths
#
# Expected structure:
#
# data/rsem_counts/HydraAHL/ID.genes.results
# data/rsem_counts/HydraRecolonization/ID.genes.results
# ------------------------------------------------------------

all_files <- file.path(
  input_dir,
  all_samples$EXP,
  paste0(
    all_samples$ID,
    ".genes.results"
  )
)

names(all_files) <- all_samples$ID


# ------------------------------------------------------------
# Check all files exist
# ------------------------------------------------------------

missing_files <- all_files[
  !file.exists(all_files)
]

if(length(missing_files) > 0) {
  
  stop(
    "Missing RSEM files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


cat(
  "[*] Importing ",
  length(all_files),
  " samples together\n",
  sep = ""
)


# ============================================================
# Import ALL samples together
# ============================================================

txi_all <- tximport(
  all_files,
  type = "rsem",
  txIn = FALSE,
  txOut = FALSE
)


# ============================================================
# Combined metadata
# ============================================================

coldata_all <- data.frame(
  
  experiment = factor(
    all_samples$EXP
  ),
  
  condition = factor(
    all_samples$treatment
  ),
  
  batch = factor(
    all_samples$batch
  ),
  
  sample_name = all_samples$sample,
  
  row.names = all_samples$ID
)


# ------------------------------------------------------------
# Match metadata order to count matrix
# ------------------------------------------------------------

coldata_all <- coldata_all[
  colnames(txi_all$counts),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(coldata_all),
    colnames(txi_all$counts)
  )
)


# ============================================================
# Create one combined DESeq2 object
#
# This object is ONLY for transformation / visualization.
# Differential expression is still done separately above.
# ============================================================

dds_combined_samples <- DESeqDataSetFromTximport(
  txi_all,
  colData = coldata_all,
  design = ~ 1
)


# ============================================================
# Count filter across combined dataset
# ============================================================

keep <- rowSums(
  counts(dds_combined_samples)
) >= count_threshold


cat(
  "[*] Genes before filtering: ",
  nrow(dds_combined_samples),
  "\n",
  sep = ""
)


dds_combined_samples <- dds_combined_samples[
  keep,
]


cat(
  "[*] Genes after filtering: ",
  nrow(dds_combined_samples),
  "\n",
  sep = ""
)


# ============================================================
# ONE common variance-stabilizing transformation
#
# This is important:
#
# HydraAHL and HydraRecolonization are transformed TOGETHER,
# putting all individual samples into the same expression space.
# ============================================================

vsd_combined_samples <- vst(
  dds_combined_samples,
  blind = TRUE
)


expr_combined <- assay(
  vsd_combined_samples
)


cat(
  "[*] Combined expression matrix: ",
  nrow(expr_combined),
  " genes x ",
  ncol(expr_combined),
  " samples\n",
  sep = ""
)


# ============================================================
# Sample-to-sample Pearson correlations
# ============================================================

sample_cor <- cor(
  expr_combined,
  method = "pearson"
)


# Avoid numerical values slightly outside [-1, 1]
sample_cor <- pmin(
  1,
  pmax(
    -1,
    sample_cor
  )
)


# Save correlation matrix
write.csv(
  sample_cor,
  file = file.path(
    combined_dir,
    "combined_samples_Pearson_correlations.csv"
  )
)


# ============================================================
# Convert correlation to distance
#
# identical profiles:
# correlation = 1
# distance    = 0
#
# unrelated:
# correlation = 0
# distance    = 1
#
# opposite:
# correlation = -1
# distance    = 2
# ============================================================

sample_distance <- as.dist(
  1 - sample_cor
)


# ============================================================
# ONE PCoA using ALL individual samples
# ============================================================

sample_pcoa <- cmdscale(
  sample_distance,
  k = 2,
  eig = TRUE,
  add = TRUE
)


# ============================================================
# Create coordinate dataframe
# ============================================================

sample_pcoa_df <- data.frame(
  
  ID = rownames(
    sample_pcoa$points
  ),
  
  PCoA1 = sample_pcoa$points[, 1],
  
  PCoA2 = sample_pcoa$points[, 2],
  
  stringsAsFactors = FALSE
)


# ============================================================
# Add metadata
# ============================================================

sample_pcoa_df$experiment <- coldata_all[
  sample_pcoa_df$ID,
  "experiment"
]

sample_pcoa_df$condition <- coldata_all[
  sample_pcoa_df$ID,
  "condition"
]

sample_pcoa_df$batch <- coldata_all[
  sample_pcoa_df$ID,
  "batch"
]

sample_pcoa_df$sample_name <- coldata_all[
  sample_pcoa_df$ID,
  "sample_name"
]


# ============================================================
# Percentage represented by PCoA axes
# ============================================================

positive_eigenvalues <- sample_pcoa$eig[
  sample_pcoa$eig > 0
]


sample_pcoa_variance <- (
  100 *
    sample_pcoa$eig /
    sum(positive_eigenvalues)
)


# ============================================================
# Save PCoA coordinates
# ============================================================

write.csv(
  sample_pcoa_df,
  file = file.path(
    combined_dir,
    "combined_sample_PCoA_coordinates.csv"
  ),
  row.names = FALSE
)


# ============================================================
# Combined PCoA plot
#
# COLOR = treatment / biological condition
# SHAPE = experiment
# ============================================================

p_sample_pcoa <- ggplot(
  sample_pcoa_df,
  
  aes(
    x = PCoA1,
    y = PCoA2,
    color = condition,
    shape = experiment
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(
      label = sample_name
    ),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  
  labs(
    
    title = "Combined sample-level PCoA",
    
    subtitle = paste0(
      "HydraAHL + HydraRecolonization | ",
      "distance = 1 - Pearson correlation"
    ),
    
    x = paste0(
      "PCoA1 (",
      round(
        sample_pcoa_variance[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PCoA2 (",
      round(
        sample_pcoa_variance[2],
        1
      ),
      "%)"
    ),
    
    color = "Treatment",
    
    shape = "Experiment"
  ) +
  
  theme_classic()


print(
  p_sample_pcoa
)


# ============================================================
# Save plot
# ============================================================

ggsave(
  
  filename = file.path(
    combined_dir,
    "HydraAHL_HydraRecolonization_combined_sample_PCoA.png"
  ),
  
  plot = p_sample_pcoa,
  
  width = 11,
  height = 8,
  dpi = 200
)


cat(
  "\n[*] Combined sample PCoA finished\n"
)

cat(
  "[*] Output: ",
  file.path(
    combined_dir,
    "HydraAHL_HydraRecolonization_combined_sample_PCoA.png"
  ),
  "\n",
  sep = ""
)

# ============================================================
# Combined contrast-level analysis
# ============================================================

cat("\n")
cat("============================================================\n")
cat("COMBINED CONTRAST-LEVEL ANALYSIS\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Validate number of contrasts
# ------------------------------------------------------------

if (length(contrast_results) < 3) {
  stop(
    "Only ",
    length(contrast_results),
    " valid contrasts were generated. ",
    "At least 3 are required for a useful 2D PCA/PCoA."
  )
}

cat(
  "[*] Contrasts available: ",
  length(contrast_results),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Genes present in every contrast
# ------------------------------------------------------------

common_genes <- Reduce(
  intersect,
  lapply(
    contrast_results,
    function(x) x$gene
  )
)

if (length(common_genes) == 0) {
  stop(
    "No genes are shared by all DESeq2 contrast tables."
  )
}

cat(
  "[*] Genes present in every contrast: ",
  length(common_genes),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Build genes x contrasts log2FC matrix
# ------------------------------------------------------------

lfc_list <- lapply(
  contrast_results,
  function(x) {

    idx <- match(
      common_genes,
      x$gene
    )

    values <- x$log2FoldChange[idx]
    names(values) <- common_genes

    values
  }
)

lfc_mat <- do.call(
  cbind,
  lfc_list
)

rownames(lfc_mat) <- common_genes
colnames(lfc_mat) <- names(contrast_results)

storage.mode(lfc_mat) <- "double"

cat(
  "[*] Initial matrix: ",
  nrow(lfc_mat),
  " genes x ",
  ncol(lfc_mat),
  " contrasts\n",
  sep = ""
)


# ------------------------------------------------------------
# Remove genes with missing/infinite LFC in any contrast
# ------------------------------------------------------------

finite_rows <- apply(
  lfc_mat,
  1,
  function(x) all(is.finite(x))
)

lfc_mat <- lfc_mat[
  finite_rows,
  ,
  drop = FALSE
]

cat(
  "[*] Genes after finite-LFC filter: ",
  nrow(lfc_mat),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Optional effect-size filter
# ------------------------------------------------------------

if (combined_min_abs_log2fc > 0) {

  keep_response <- apply(
    abs(lfc_mat),
    1,
    max
  ) >= combined_min_abs_log2fc

  lfc_mat <- lfc_mat[
    keep_response,
    ,
    drop = FALSE
  ]

  cat(
    "[*] Genes after |log2FC| >= ",
    combined_min_abs_log2fc,
    " in at least one contrast: ",
    nrow(lfc_mat),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# Remove genes with zero variance between contrasts
# ------------------------------------------------------------

gene_sd <- apply(
  lfc_mat,
  1,
  sd
)

variable_rows <- is.finite(gene_sd) & gene_sd > 0

lfc_mat <- lfc_mat[
  variable_rows,
  ,
  drop = FALSE
]

if (nrow(lfc_mat) < 2) {
  stop(
    "Too few variable genes remain for combined analysis."
  )
}

cat(
  "[*] Genes used in combined analysis: ",
  nrow(lfc_mat),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Save combined LFC matrix
# ------------------------------------------------------------

write.csv(
  data.frame(
    gene = rownames(lfc_mat),
    lfc_mat,
    check.names = FALSE
  ),
  file = file.path(
    combined_dir,
    "combined_log2FC_matrix.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Combined contrast metadata
# ------------------------------------------------------------

metadata <- do.call(
  rbind,
  contrast_metadata
)

rownames(metadata) <- NULL

metadata <- metadata[
  match(
    colnames(lfc_mat),
    metadata$contrast
  ),
  ,
  drop = FALSE
]

if (any(is.na(metadata$contrast))) {
  stop(
    "Could not match combined-analysis columns to contrast metadata."
  )
}

write.csv(
  metadata,
  file = file.path(
    combined_dir,
    "contrast_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PCA of DESeq2 contrasts
# ============================================================

cat("[*] Running PCA of contrast log2FC profiles\n")

pca <- prcomp(
  t(lfc_mat),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  contrast = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_df <- merge(
  pca_df,
  metadata,
  by = "contrast",
  all.x = TRUE,
  sort = FALSE
)

pca_df <- pca_df[
  match(
    rownames(pca$x),
    pca_df$contrast
  ),
  ,
  drop = FALSE
]

write.csv(
  pca_df,
  file = file.path(
    combined_dir,
    "combined_PCA_coordinates.csv"
  ),
  row.names = FALSE
)

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = experiment,
    label = contrast
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of transcriptional responses",
    subtitle = paste0(
      "DESeq2 log2FC profiles; ",
      nrow(lfc_mat),
      " genes"
    ),
    x = paste0(
      "PC1 (",
      round(pca_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2], 1),
      "%)"
    ),
    color = "Experiment"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    combined_dir,
    "combined_log2FC_PCA.png"
  ),
  plot = p_pca,
  width = 10,
  height = 7,
  dpi = 200
)


# ------------------------------------------------------------
# PCA gene loadings
# ------------------------------------------------------------

loadings <- data.frame(
  gene = rownames(pca$rotation),
  PC1 = pca$rotation[, 1],
  PC2 = pca$rotation[, 2],
  stringsAsFactors = FALSE
)

write.csv(
  loadings,
  file = file.path(
    combined_dir,
    "PCA_gene_loadings.csv"
  ),
  row.names = FALSE
)

top_pc1 <- loadings[
  order(
    abs(loadings$PC1),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

top_pc2 <- loadings[
  order(
    abs(loadings$PC2),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  head(top_pc1, 50),
  file = file.path(
    combined_dir,
    "top50_PC1_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  head(top_pc2, 50),
  file = file.path(
    combined_dir,
    "top50_PC2_genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# Pearson correlation between contrasts
# ============================================================

cat("[*] Calculating Pearson correlations between contrasts\n")

cor_mat <- cor(
  lfc_mat,
  method = "pearson"
)

# Guard against tiny numerical values outside [-1, 1].
cor_mat <- pmin(
  1,
  pmax(-1, cor_mat)
)

write.csv(
  cor_mat,
  file = file.path(
    combined_dir,
    "contrast_Pearson_correlations.csv"
  )
)


# ------------------------------------------------------------
# Correlation-based clustering
# ------------------------------------------------------------

cor_distance <- as.dist(
  1 - cor_mat
)

hc <- hclust(
  cor_distance,
  method = "average"
)

contrast_order <- hc$labels[
  hc$order
]

cor_ordered <- cor_mat[
  contrast_order,
  contrast_order,
  drop = FALSE
]


# ------------------------------------------------------------
# Correlation heatmap with ggplot2
# ------------------------------------------------------------

cor_long <- as.data.frame(
  as.table(cor_ordered),
  stringsAsFactors = FALSE
)

names(cor_long) <- c(
  "Contrast1",
  "Contrast2",
  "Correlation"
)

cor_long$Contrast1 <- factor(
  cor_long$Contrast1,
  levels = contrast_order
)

cor_long$Contrast2 <- factor(
  cor_long$Contrast2,
  levels = rev(contrast_order)
)

p_cor <- ggplot(
  cor_long,
  aes(
    x = Contrast1,
    y = Contrast2,
    fill = Correlation
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson\ncorrelation"
  ) +
  labs(
    title = "Similarity of transcriptional responses",
    subtitle = "Pearson correlation of DESeq2 log2 fold changes",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    combined_dir,
    "contrast_correlation_heatmap.png"
  ),
  plot = p_cor,
  width = 10,
  height = 8,
  dpi = 200
)


# ============================================================
# PCoA using distance = 1 - Pearson correlation
# ============================================================

cat("[*] Running correlation-based PCoA\n")

pcoa <- cmdscale(
  cor_distance,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_df <- data.frame(
  contrast = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)

pcoa_df <- merge(
  pcoa_df,
  metadata,
  by = "contrast",
  all.x = TRUE,
  sort = FALSE
)

pcoa_df <- pcoa_df[
  match(
    rownames(pcoa$points),
    pcoa_df$contrast
  ),
  ,
  drop = FALSE
]

positive_eigenvalues <- pcoa$eig[
  pcoa$eig > 0
]

if (length(positive_eigenvalues) > 0) {
  pcoa_variance <- 100 * pcoa$eig / sum(positive_eigenvalues)
} else {
  pcoa_variance <- rep(NA_real_, length(pcoa$eig))
}

write.csv(
  pcoa_df,
  file = file.path(
    combined_dir,
    "combined_PCoA_coordinates.csv"
  ),
  row.names = FALSE
)

p_pcoa <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = experiment,
    label = contrast
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "PCoA of transcriptional responses",
    subtitle = "Distance = 1 - Pearson correlation of DESeq2 log2FC profiles",
    x = paste0(
      "PCoA1 (",
      round(pcoa_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(pcoa_variance[2], 1),
      "%)"
    ),
    color = "Experiment"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    combined_dir,
    "combined_log2FC_PCoA.png"
  ),
  plot = p_pcoa,
  width = 10,
  height = 7,
  dpi = 200
)


# ============================================================
# Final summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat(
  "Experiments successfully analysed: ",
  length(dds_objects),
  "\n",
  sep = ""
)
cat(
  "DESeq2 contrasts successfully analysed: ",
  length(contrast_results),
  "\n",
  sep = ""
)
cat(
  "Genes in combined analysis: ",
  nrow(lfc_mat),
  "\n",
  sep = ""
)
cat(
  "Combined outputs: ",
  combined_dir,
  "\n",
  sep = ""
)
cat("============================================================\n")
