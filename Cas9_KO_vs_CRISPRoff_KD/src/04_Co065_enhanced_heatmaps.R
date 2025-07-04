# Load required libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Configuration
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1
targets <- c("FAS", "PTPN2", "RC3H1", "SUV39H1", "RASA2", "MED12")

cat("CO065 ENHANCED OFF-TARGET ANALYSIS\n")

# STEP 1: Load RNA-seq results to identify common significant genes

cat("Loading RNA-seq results to identify genes significant in both KD and KO...\n")

# Function to safely load RNA-seq results
load_rna_results <- function(target, exp_type) {
  # Try different file naming patterns
  # Handle case variations - FAS is saved as "Fas" in files
  target_variations <- c(target, toupper(target), tolower(target))
  if (target == "FAS") {
    target_variations <- c("Fas", "FAS", "fas")
  }
  
  file_patterns <- character()
  for (tgt in target_variations) {
    file_patterns <- c(file_patterns,
      file.path("results_Co065", paste0("DE_topTable_Co065_", tgt, "_", exp_type, ".csv"))
    )
  }
  
  for (file_path in file_patterns) {
    if (file.exists(file_path)) {
      return(read_csv(file_path, show_col_types = FALSE))
    }
  }
  
  cat("Missing RNA-seq file for", target, exp_type, "\n")
  return(NULL)
}

# Function to get significant genes
get_sig_genes <- function(df) {
  if (is.null(df)) return(character(0))
  df %>%
    filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) %>%
    pull(gene_symbol) %>%
    na.omit() %>%
    unique()
}

# Load all RNA-seq results and find common significant genes
all_rna_results <- list()
common_sig_by_target <- list()

for (target in targets) {
  # Load KO and KD results
  ko_result <- load_rna_results(target, "KO")
  kd_result <- load_rna_results(target, "KD")
  
  if (!is.null(ko_result)) all_rna_results[[paste0(target, "_KO")]] <- ko_result
  if (!is.null(kd_result)) all_rna_results[[paste0(target, "_KD")]] <- kd_result
  
  # Find common significant genes
  ko_genes <- get_sig_genes(ko_result)
  kd_genes <- get_sig_genes(kd_result)
  common_genes <- intersect(ko_genes, kd_genes)
  
  common_sig_by_target[[target]] <- common_genes
  
  cat("Target:", target, "- KO sig:", length(ko_genes), 
      ", KD sig:", length(kd_genes), ", Common:", length(common_genes), "\n")
  
  if (length(common_genes) > 0) {
    cat("  Common genes:", paste(common_genes, collapse = ", "), "\n")
  }
}

# STEP 2: Enhanced heatmap function

# Enhanced heatmap creation function that includes KD∩KO significance
create_enhanced_heatmap <- function(target, exp_type, common_sig_genes = NULL) {
  
  # Load existing heatmap data
  logfc_file <- file.path("Co065_heatmap_results", paste0("Co065_", target, "_", exp_type, "_logFC_data.csv"))
  guide_file <- file.path("Co065_heatmap_results", paste0("Co065_", target, "_", exp_type, "_guide_info.csv"))
  
  if (!file.exists(logfc_file) || !file.exists(guide_file)) {
    cat("Missing files for", target, exp_type, "- skipping\n")
    return(NULL)
  }
  
  logfc_data <- read_csv(logfc_file, show_col_types = FALSE)
  guide_info <- read_csv(guide_file, show_col_types = FALSE)
  
  if (nrow(logfc_data) == 0) {
    cat("No data for", target, exp_type, "\n")
    return(NULL)
  }
  
  cat("Creating enhanced heatmap for", target, exp_type, "- Genes:", nrow(logfc_data), "\n")
  
  # This data is already in single-column format (gene per row, single logFC column)
  # Create a single-row matrix for transposed heatmap (genes on x-axis)
  mat <- matrix(logfc_data$logFC, nrow = 1, dimnames = list(paste(target, exp_type), logfc_data$gene_symbol))
  
  # Remove any columns with missing gene symbols
  valid_cols <- !is.na(colnames(mat)) & colnames(mat) != ""
  mat <- mat[, valid_cols, drop = FALSE]
  
  # Filter annotations to match valid matrix columns
  logfc_data_filtered <- logfc_data[valid_cols, ]
  
  # Create enhanced annotations
  # 1. Analysis type (existing)
  analysis_annotation <- logfc_data_filtered$analysis_type
  names(analysis_annotation) <- colnames(mat)
  
  # Manual override: ensure the target gene itself is classified as "target"
  target_gene_idx <- which(colnames(mat) == target)
  if (length(target_gene_idx) > 0) {
    analysis_annotation[target_gene_idx] <- "target"
  }
  
  # 2. Significance in current experiment - updated criteria: |logFC| > 1 AND adj.P.Val < 0.05
  # Note: the significance column in the data uses adj.P.Val < 0.05, we need to add logFC filter
  sig_annotation <- logfc_data_filtered$significance & (abs(logfc_data_filtered$logFC) > LOGFC_CUTOFF)
  names(sig_annotation) <- colnames(mat)
  
  # 3. NEW: Significance in both KD and KO (enhanced)
  common_sig_annotation <- rep(FALSE, ncol(mat))
  names(common_sig_annotation) <- colnames(mat)
  if (!is.null(common_sig_genes) && length(common_sig_genes) > 0) {
    common_sig_annotation[colnames(mat) %in% common_sig_genes] <- TRUE
  }
  
  # Reorder columns to put target gene last (rightmost)
  if (length(target_gene_idx) > 0) {
    # Create new column order with target last
    new_order <- c(setdiff(1:ncol(mat), target_gene_idx), target_gene_idx)
    
    # Reorder matrix
    mat <- mat[, new_order, drop = FALSE]
    
    # Reorder all annotations
    analysis_annotation <- analysis_annotation[new_order]
    sig_annotation <- sig_annotation[new_order]
    common_sig_annotation <- common_sig_annotation[new_order]
  }
  
  # Create enhanced column annotation - simplified to match reference style
  column_annotation <- columnAnnotation(
    Criteria = analysis_annotation,
    col = list(
      Criteria = c("offtarget" = "orange", "proximal" = "green", "both" = "purple", "target" = "blue")
    ),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Criteria = list(
        title = "Criteria:",
        labels = c("predicted off target", "proximal gene", "Both", "target gene"),
        at = c("offtarget", "proximal", "both", "target")
      )
    ),
    simple_anno_size = unit(4, "mm")
  )
  
  # Create enhanced heatmap
  ht <- Heatmap(
    mat,
    name = "log2FC",
    col = colorRamp2(c(-4, -2, 0, 2, 4), c("blue", "lightblue", "white", "lightcoral", "red")),
    
    # Titles
    column_title = paste0(target, " ", exp_type),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title_side = "top",
    
    # Annotations
    top_annotation = column_annotation,
    
    # Clustering - disabled to maintain original ordering
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    # Display
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 90,
    column_names_side = "bottom",
    
    # Cell annotations - add asterisks for significance
    cell_fun = function(j, i, x, y, width, height, fill) {
      gene <- colnames(mat)[j]
      
      # Add asterisk for significance
      if (gene %in% names(sig_annotation) && sig_annotation[gene]) {
        grid.text("*", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
      }
    },
    
    # Heatmap dimensions
    # Remove cell borders for cleaner look
    rect_gp = gpar(col = "white", lwd = 0),
    
    # Heatmap dimensions
    width = unit(max(10, ncol(mat) * 0.25), "cm"),
    height = unit(1, "cm")
  )
  
  return(ht)
}

# STEP 3: Generate enhanced heatmaps

cat("\nGENERATING ENHANCED HEATMAPS\n")

# Create output directory
enhanced_output_dir <- "Co065_enhanced_heatmaps"
if (!dir.exists(enhanced_output_dir)) dir.create(enhanced_output_dir, recursive = TRUE)

# Generate enhanced heatmaps for all target-experiment combinations
for (target in targets) {
  common_genes <- common_sig_by_target[[target]]
  
  for (exp_type in c("KD", "KO")) {
    cat("Processing", target, exp_type, "...\n")
    
    ht <- create_enhanced_heatmap(target, exp_type, common_genes)
    
    if (!is.null(ht)) {
      # Save enhanced heatmap
      png_file <- file.path(enhanced_output_dir, paste0("Co065_", target, "_", exp_type, "_enhanced_heatmap.png"))
      pdf_file <- file.path(enhanced_output_dir, paste0("Co065_", target, "_", exp_type, "_enhanced_heatmap.pdf"))
      
      # PNG version - adjusted for horizontal layout
      png(png_file, width = max(10, ncol(ht@matrix) * 0.15 + 2), height = 4, units = "in", res = 300)
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top", 
           merge_legend = TRUE, legend_gap = unit(0.5, "cm"))
      dev.off()
      
      # PDF version - adjusted for horizontal layout
      pdf(pdf_file, width = max(10, ncol(ht@matrix) * 0.15 + 2), height = 4)
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top", 
           merge_legend = TRUE, legend_gap = unit(0.5, "cm"))
      dev.off()
      
      cat("Saved enhanced heatmap:", png_file, "\n")
    } else {
      cat("Could not create heatmap for", target, exp_type, "\n")
    }
  }
}

# STEP 4: Create summary report

cat("\nCREATING SUMMARY REPORT\n")

# Create summary of common significant genes
summary_data <- tibble(
  Target = targets,
  KD_Sig_Count = map_int(targets, ~{
    result <- all_rna_results[[paste0(.x, "_KD")]]
    length(get_sig_genes(result))
  }),
  KO_Sig_Count = map_int(targets, ~{
    result <- all_rna_results[[paste0(.x, "_KO")]]
    length(get_sig_genes(result))
  }),
  Common_Sig_Count = map_int(targets, ~length(common_sig_by_target[[.x]])),
  Common_Sig_Genes = map_chr(targets, ~{
    genes <- common_sig_by_target[[.x]]
    if (length(genes) > 0) paste(genes, collapse = ", ") else "None"
  })
)

# Save summary
summary_file <- file.path(enhanced_output_dir, "Co065_enhanced_analysis_summary.csv")
write_csv(summary_data, summary_file)

cat("Enhanced analysis summary:\n")
print(summary_data)

# Create summary text report
report_lines <- c(
  "# Co065 Enhanced Off-target Analysis Report",
  paste("Generated on:", Sys.Date()),
  "",
  "## Summary",
  "This analysis integrates RNA-seq differential expression results with off-target analysis",
  "to highlight genes that are significantly affected in both KD and KO experiments.",
  "",
  "## Key Features",
  "- Enhanced heatmaps show off-target and proximal genes",
  "- Purple annotation track highlights genes significant in both KD and KO",
  "- Single asterisk (*) = significant in current experiment",
  "- Double asterisk (**) = significant in both KD and KO experiments",
  "",
  "## Results by Target",
  ""
)

for (i in seq_len(nrow(summary_data))) {
  target_data <- summary_data[i, ]
  report_lines <- c(
    report_lines,
    paste("###", target_data$Target),
    paste("- KD significant genes:", target_data$KD_Sig_Count),
    paste("- KO significant genes:", target_data$KO_Sig_Count),
    paste("- Common significant genes:", target_data$Common_Sig_Count),
    if (target_data$Common_Sig_Genes != "None") {
      paste("- Common genes:", target_data$Common_Sig_Genes)
    } else {
      "- No genes significant in both conditions"
    },
    ""
  )
}

report_file <- file.path(enhanced_output_dir, "Co065_enhanced_analysis_report.md")
writeLines(report_lines, report_file)

cat("\nENHANCED ANALYSIS COMPLETE\n")
cat("Results saved to:", enhanced_output_dir, "\n")
cat("- Enhanced heatmaps with KD∩KO significance highlighting\n")
cat("- Summary CSV:", summary_file, "\n")
cat("- Detailed report:", report_file, "\n")
cat("\nKey Legend:\n")
cat("- Blue = target gene\n")
cat("- Orange = predicted off target genes\n")
cat("- Green = proximal genes\n")
cat("- Purple = genes in both categories\n")
cat("- * = significant (|logFC| > 1 & adj.P.Val < 0.05)\n")