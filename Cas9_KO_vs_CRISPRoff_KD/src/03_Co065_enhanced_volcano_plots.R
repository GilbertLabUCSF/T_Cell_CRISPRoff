library(tidyverse)
library(ggrepel)
library(patchwork)

# Configuration
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1
targets <- c("FAS", "PTPN2", "RC3H1", "SUV39H1", "RASA2", "MED12")

# Function to safely load RNA-seq results
load_rna_results <- function(target, exp_type) {
  file_patterns <- c(
    file.path("results_Co065", paste0("DE_topTable_Co065_", target, "_", exp_type, ".csv")),
    file.path("results_Co065", paste0("DE_topTable_Co065_", toupper(target), "_", exp_type, ".csv"))
  )
  
  for (file_path in file_patterns) {
    if (file.exists(file_path)) {
      return(read_csv(file_path, show_col_types = FALSE))
    }
  }
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
  ko_result <- load_rna_results(target, "KO")
  kd_result <- load_rna_results(target, "KD")
  
  if (!is.null(ko_result)) all_rna_results[[paste0(target, "_KO")]] <- ko_result
  if (!is.null(kd_result)) all_rna_results[[paste0(target, "_KD")]] <- kd_result
  
  ko_genes <- get_sig_genes(ko_result)
  kd_genes <- get_sig_genes(kd_result)
  common_genes <- intersect(ko_genes, kd_genes)
  
  common_sig_by_target[[target]] <- common_genes
  
  cat("Target:", target, "- Common significant genes:", length(common_genes), "\n")
}

create_enhanced_volcano_plot <- function(df_results, contrast_id, target_gene,
                                       common_genes = NULL,
                                       logfc_cutoff = LOGFC_CUTOFF,
                                       p_cutoff = P_VALUE_CUTOFF) {
  
  if (is.null(df_results) || nrow(df_results) == 0) return(NULL)
  
  # Create significance categories with special highlighting for common genes
  df_results <- df_results %>%
    mutate(Significance = case_when(
      !is.na(gene_symbol) & gene_symbol == target_gene ~ "Target Gene",
      !is.null(common_genes) & gene_symbol %in% common_genes & logFC > 0 ~ "Up-regulated",
      !is.null(common_genes) & gene_symbol %in% common_genes & logFC < 0 ~ "Down-regulated",
      adj.P.Val < p_cutoff & logFC >  logfc_cutoff ~ "Up-regulated",
      adj.P.Val < p_cutoff & logFC < -logfc_cutoff ~ "Down-regulated",
      TRUE ~ "Not Significant"
    )) %>%
    mutate(
      # Track which genes are shared between KD and KO for labeling
      is_shared = !is.null(common_genes) & gene_symbol %in% common_genes,
      Significance = factor(Significance, 
                           levels = c("Target Gene", "Up-regulated", "Down-regulated", "Not Significant"))
    )
  
  # Create the plot with proper layering to ensure target gene is on top
  p <- ggplot(df_results, aes(logFC, -log10(adj.P.Val))) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", colour = "red", alpha = 0.7) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
               linetype = "dashed", colour = "grey50") +
    
    # Plot all points with solid colors (no transparency)
    geom_point(aes(colour = Significance, 
                   size = ifelse(Significance == "Target Gene", 2.5, 1.5))) +
    
    # Plot target gene again on top to ensure maximum visibility
    geom_point(data = df_results %>% filter(Significance == "Target Gene"),
               aes(colour = Significance), size = 3) +
    
    scale_colour_manual(
      values = c(
        "Target Gene" = "orange",
        "Up-regulated" = "red", 
        "Down-regulated" = "blue", 
        "Not Significant" = "grey60"
      ),
      labels = c(
        "Target Gene" = "Target Gene",
        "Up-regulated" = "Up-regulated", 
        "Down-regulated" = "Down-regulated", 
        "Not Significant" = "Not Significant"
      )
    ) +
    scale_size_identity() +
    
    # Label target gene
    ggrepel::geom_text_repel(
      data = df_results %>% filter(Significance == "Target Gene"),
      aes(label = gene_symbol),
      box.padding = 0.5, point.padding = 0.3,
      max.overlaps = Inf, size = 3, fontface = "bold",
      color = "orange"
    ) +
    
    # Label genes significant in both KD and KO (shared genes)
    ggrepel::geom_text_repel(
      data = df_results %>% filter(is_shared == TRUE),
      aes(label = gene_symbol, color = Significance),
      box.padding = 0.4, point.padding = 0.15,
      max.overlaps = 15, size = 2.5, fontface = "bold",
      show.legend = FALSE
    ) +
    
    labs(title = paste(contrast_id, "- Enhanced Volcano Plot"),
         subtitle = paste("Labeled genes are significant in both KD & KO (n =", 
                         length(common_genes %||% 0), ")"),
         x = "log2 Fold Change",
         y = "-log10 Adjusted P-value") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  return(p)
}


# Create output directory
enhanced_output_dir <- "Co065_enhanced_volcano_plots"
if (!dir.exists(enhanced_output_dir)) dir.create(enhanced_output_dir, recursive = TRUE)

# Generate enhanced volcano plots for each target
all_plots <- list()

for (target in targets) {
  common_genes <- common_sig_by_target[[target]]
  
  for (exp_type in c("KD", "KO")) {
    contrast_name <- paste0(target, "_", exp_type)
    
    if (contrast_name %in% names(all_rna_results)) {
      cat("Creating enhanced volcano plot for", contrast_name, "...\n")
      
      volcano_plot <- create_enhanced_volcano_plot(
        all_rna_results[[contrast_name]], 
        contrast_name, 
        target,
        common_genes
      )
      
      if (!is.null(volcano_plot)) {
        # Save individual plot
        png_file <- file.path(enhanced_output_dir, paste0("Co065_", contrast_name, "_enhanced_volcano.png"))
        ggsave(png_file, volcano_plot, width = 8, height = 6, dpi = 300)
        cat("Saved:", png_file, "\n")
        
        # Store for combined plot
        all_plots[[contrast_name]] <- volcano_plot
      }
    }
  }
}


for (target in targets) {
  kd_plot_name <- paste0(target, "_KD")
  ko_plot_name <- paste0(target, "_KO")
  
  if (kd_plot_name %in% names(all_plots) && ko_plot_name %in% names(all_plots)) {
    # Get the y-axis range from both plots to set shared limits
    kd_data <- all_plots[[kd_plot_name]]$data
    ko_data <- all_plots[[ko_plot_name]]$data
    
    # Calculate shared y-axis limits
    y_max <- max(c(-log10(kd_data$adj.P.Val), -log10(ko_data$adj.P.Val)), na.rm = TRUE)
    y_min <- 0
    
    # Update both plots with shared y-axis limits
    kd_plot_shared <- all_plots[[kd_plot_name]] + 
      ylim(y_min, y_max * 1.05)  # Add 5% padding
    
    ko_plot_shared <- all_plots[[ko_plot_name]] + 
      ylim(y_min, y_max * 1.05)  # Add 5% padding
    
    combined_plot <- kd_plot_shared + ko_plot_shared +
      plot_layout(ncol = 2, guides = "collect")
    
    combined_file <- file.path(enhanced_output_dir, paste0("Co065_", target, "_combined_enhanced_volcano.png"))
    ggsave(combined_file, combined_plot, width = 14, height = 6, dpi = 300)
    cat("Saved combined plot:", combined_file, "\n")
  }
}


# Filter plots that exist
existing_plots <- all_plots[names(all_plots) %in% names(all_rna_results)]

if (length(existing_plots) > 0) {
  # Create mega combined plot
  n_plots <- length(existing_plots)
  ncols <- 3
  nrows <- ceiling(n_plots / ncols)
  
  mega_plot <- wrap_plots(existing_plots, ncol = ncols, guides = "collect")
  
  mega_file <- file.path(enhanced_output_dir, "Co065_all_targets_enhanced_volcano_plots.png")
  ggsave(mega_file, mega_plot, 
         width = 5 * ncols, height = 4 * nrows, dpi = 300)
  cat("Saved mega plot:", mega_file, "\n")
}

# Create summary of common genes per target
summary_stats <- tibble(
  Target = targets,
  KD_Available = map_lgl(targets, ~paste0(.x, "_KD") %in% names(all_rna_results)),
  KO_Available = map_lgl(targets, ~paste0(.x, "_KO") %in% names(all_rna_results)),
  Common_Genes_Count = map_int(targets, ~length(common_sig_by_target[[.x]])),
  Top_Common_Genes = map_chr(targets, ~{
    genes <- common_sig_by_target[[.x]]
    if (length(genes) > 0) {
      paste(head(genes, 5), collapse = ", ")
    } else {
      "None"
    }
  })
)

# Save summary
summary_file <- file.path(enhanced_output_dir, "Co065_enhanced_volcano_summary.csv")
write_csv(summary_stats, summary_file)