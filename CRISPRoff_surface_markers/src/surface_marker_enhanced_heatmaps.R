# Surface Marker Enhanced Heatmap Generation
# Adapted from Co062 enhanced heatmaps
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Configuration
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1
targets <- c("CD55", "CD81")

cat("SURFACE MARKER ENHANCED HEATMAP GENERATION\n")
cat("==========================================\n\n")

# Enhanced heatmap creation function (adapted from Co062)
create_enhanced_heatmap <- function(guide_name) {
  
  # Extract target from guide name
  target <- strsplit(guide_name, "_")[[1]][1]
  
  # Load existing heatmap data
  logfc_file <- file.path("heatmap_results", paste0(guide_name, "_logFC_data.csv"))
  guide_file <- file.path("heatmap_results", paste0(guide_name, "_guide_info.csv"))
  
  if (!file.exists(logfc_file) || !file.exists(guide_file)) {
    cat("Missing files for", guide_name, "- skipping\n")
    return(NULL)
  }
  
  logfc_data <- read_csv(logfc_file, show_col_types = FALSE)
  guide_info <- read_csv(guide_file, show_col_types = FALSE)
  
  if (nrow(logfc_data) == 0) {
    cat("No data for", guide_name, "\n")
    return(NULL)
  }
  
  cat("Creating enhanced heatmap for", guide_name, "- Genes:", nrow(logfc_data), "\n")
  
  # Create single-row matrix for transposed heatmap (genes on x-axis)
  mat <- matrix(logfc_data$logFC, nrow = 1, dimnames = list(guide_name, logfc_data$gene_symbol))
  
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
  sig_annotation <- logfc_data_filtered$significance & (abs(logfc_data_filtered$logFC) > LOGFC_CUTOFF)
  names(sig_annotation) <- colnames(mat)
  
  # Reorder columns to put target gene last (rightmost)
  if (length(target_gene_idx) > 0) {
    # Create new column order with target last
    new_order <- c(setdiff(1:ncol(mat), target_gene_idx), target_gene_idx)
    
    # Reorder matrix
    mat <- mat[, new_order, drop = FALSE]
    
    # Reorder all annotations
    analysis_annotation <- analysis_annotation[new_order]
    sig_annotation <- sig_annotation[new_order]
  }
  
  # Create enhanced column annotation
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
    col = colorRamp2(c(-8, -4, 0, 4, 8), c("blue", "lightblue", "white", "lightcoral", "red")),
    
    # Titles
    column_title = paste0(target, " CRISPRoff (3 pooled guides)"),
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
    
    # Remove cell borders for cleaner look
    rect_gp = gpar(col = "white", lwd = 0),
    
    # Heatmap dimensions
    width = unit(max(10, ncol(mat) * 0.25), "cm"),
    height = unit(1, "cm")
  )
  
  return(ht)
}

# Create output directory
enhanced_output_dir <- "enhanced_heatmaps"
if (!dir.exists(enhanced_output_dir)) dir.create(enhanced_output_dir, recursive = TRUE)

# Generate enhanced heatmaps for both targets
all_guides <- paste0(targets, "_pool")

for (guide_name in all_guides) {
  cat("Processing", guide_name, "...\n")
  
  ht <- create_enhanced_heatmap(guide_name)
  
  if (!is.null(ht)) {
    # Save enhanced heatmap
    png_file <- file.path(enhanced_output_dir, paste0(guide_name, "_enhanced_heatmap.png"))
    pdf_file <- file.path(enhanced_output_dir, paste0(guide_name, "_enhanced_heatmap.pdf"))
      
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
    cat("Could not create heatmap for", guide_name, "\n")
  }
}

# Create combined heatmap
cat("\nCreating combined heatmap...\n")

all_heatmaps <- list()
for (guide_name in all_guides) {
  ht <- create_enhanced_heatmap(guide_name)
  if (!is.null(ht)) {
    all_heatmaps[[guide_name]] <- ht
  }
}

if (length(all_heatmaps) == 2) {
  combined_ht <- all_heatmaps[[1]] %v% all_heatmaps[[2]]
  
  # Save combined heatmap
  png_file <- file.path(enhanced_output_dir, "CD55_CD81_combined_enhanced_heatmap.png")
  pdf_file <- file.path(enhanced_output_dir, "CD55_CD81_combined_enhanced_heatmap.pdf")
  
  # PNG version
  png(png_file, width = max(14, max(ncol(all_heatmaps[[1]]@matrix), ncol(all_heatmaps[[2]]@matrix)) * 0.15 + 2), 
      height = 8, units = "in", res = 300)
  draw(combined_ht, heatmap_legend_side = "right", annotation_legend_side = "top", 
       merge_legend = TRUE, legend_gap = unit(0.5, "cm"))
  dev.off()
  
  # PDF version
  pdf(pdf_file, width = max(14, max(ncol(all_heatmaps[[1]]@matrix), ncol(all_heatmaps[[2]]@matrix)) * 0.15 + 2), 
      height = 8)
  draw(combined_ht, heatmap_legend_side = "right", annotation_legend_side = "top", 
       merge_legend = TRUE, legend_gap = unit(0.5, "cm"))
  dev.off()
  
  cat("Saved combined heatmap:", png_file, "\n")
}

# Create summary report
cat("\nCREATING SUMMARY REPORT\n")

# Load RNA-seq results for summary
summary_data <- tibble(
  Guide = character(),
  Target = character(),
  Sig_Count = integer()
)

for (guide_name in all_guides) {
  target <- strsplit(guide_name, "_")[[1]][1]
  
  # Load RNA-seq result
  de_file <- file.path("results_dge", paste0("DE_", target, "_vs_NTC.csv"))
  if (file.exists(de_file)) {
    de_data <- read_csv(de_file, show_col_types = FALSE)
    sig_genes <- de_data %>%
      filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) %>%
      nrow()
    
    summary_data <- summary_data %>%
      add_row(
        Guide = guide_name,
        Target = target,
        Sig_Count = sig_genes
      )
  }
}

# Save summary
summary_file <- file.path(enhanced_output_dir, "surface_marker_enhanced_analysis_summary.csv")
write_csv(summary_data, summary_file)

cat("Enhanced analysis summary:\n")
print(summary_data)

# Create summary text report
report_lines <- c(
  "# Surface Marker Enhanced Off-target Analysis Report",
  paste("Generated on:", Sys.Date()),
  "",
  "## Summary",
  "This analysis integrates RNA-seq differential expression results with off-target analysis",
  "to highlight genes that are significantly affected by each CRISPRoff guide pool.",
  "",
  "## Key Features",
  "- Enhanced heatmaps show off-target and proximal genes",
  "- Colored annotation track shows gene categories",
  "- Asterisk (*) = significant (|logFC| > 1 & adj.P.Val < 0.05)",
  "- Target genes positioned at the rightmost position",
  "",
  "## Results by Target",
  ""
)

for (i in seq_len(nrow(summary_data))) {
  target_data <- summary_data[i, ]
  report_lines <- c(
    report_lines,
    paste("###", target_data$Guide),
    paste("- Target:", target_data$Target),
    paste("- Significant genes:", target_data$Sig_Count),
    ""
  )
}

report_file <- file.path(enhanced_output_dir, "surface_marker_enhanced_analysis_report.md")
writeLines(report_lines, report_file)

cat("\nENHANCED ANALYSIS COMPLETE\n")
cat("Results saved to:", enhanced_output_dir, "\n")
cat("- Enhanced heatmaps with significance highlighting\n")
cat("- Summary CSV:", summary_file, "\n")
cat("- Detailed report:", report_file, "\n")
cat("\nKey Legend:\n")
cat("- Blue = target gene\n")
cat("- Orange = predicted off target genes\n")
cat("- Green = proximal genes\n")
cat("- Purple = genes in both categories\n")
cat("- * = significant (|logFC| > 1 & adj.P.Val < 0.05)\n")