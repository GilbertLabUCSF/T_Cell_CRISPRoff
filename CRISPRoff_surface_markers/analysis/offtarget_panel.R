library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggsignif)

# Load DE results
de_results <- read_csv("results_dge/DE_CD55_vs_NTC.csv")

# Load and prepare the normalized counts data - only CD55 and NTC samples
norm_counts_df <- read_csv("results_dge/normalized_counts.csv") %>%
  select(gene_id, gene_symbol, 
         # NTC samples
         contains("Neg"),
         # CD55 samples
         contains("55_33b")) %>%
  pivot_longer(-c(gene_id, gene_symbol), names_to = "sample", values_to = "log2_cpm") %>%
  mutate(Target = if_else(grepl("Neg", sample), "NTC", "CD55"))

# Print sample counts to verify data
cat("Number of samples by target:\n")
print(table(norm_counts_df$Target))

create_gene_expression_plot <- function(gene_symbol) {
  data_subset <- norm_counts_df %>%
    filter(gene_symbol == !!gene_symbol)
  
  if(nrow(data_subset) == 0) {
    warning(paste("No data found for gene:", gene_symbol))
    return(NULL)
  }
  
  # Get adjusted p-value for this gene
  p_val <- de_results %>%
    filter(gene_symbol == !!gene_symbol) %>%
    pull(adj.P.Val) %>%
    first()
  
  p_text <- if(!is.na(p_val)) {
    sprintf("p = %.6f", p_val)  # Show 6 decimal places
  } else "ns"
  
  ggplot(data_subset, aes(x = Target, y = log2_cpm, fill = Target)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.6) +
    scale_fill_manual(values = c("NTC" = "grey", "CD55" = "darkred")) +
    geom_signif(
      comparisons = list(c("NTC", "CD55")),
      annotations = p_text,
      y_position = 50,
      tip_length = 0
    ) +
    labs(title = paste0(gene_symbol),
         y = "log2 CPM", 
         x = "") +
    theme_pubr() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )+
    ylim(0, 55)
}

# Create plots for multiple genes
plot_gene_panel <- function(genes, ncol = 3) {
  plots <- lapply(genes, function(g) create_gene_expression_plot(g))
  valid_plots <- plots[!sapply(plots, is.null)]
  
  if(length(valid_plots) == 0) {
    stop("No valid plots were created. Check gene names.")
  }
  
  wrap_plots(valid_plots, ncol = ncol)
}

# Example usage:
genes_of_interest <- c("CD55", "EFCAB13-DT", "SMURF2")
panel <- plot_gene_panel(genes_of_interest)
ggsave("plots_dge/gene_panel.pdf", panel, width = 4*3, height = 4*ceiling(length(genes_of_interest)/3))