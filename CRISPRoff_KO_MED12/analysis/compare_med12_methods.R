library(tidyverse)
library(ggrepel)
library(ggpubr)

# Read the DE results
crispri_results <- read_csv("CRISPRoff/results_dge/DE_targetMED12_vs_NTC.csv")
ko_results <- read_csv("CRISPRoff_KO_MED12/results_med12_dge/DE_MED12_vs_AAVS1.csv")

# Combine the results
combined_results <- crispri_results %>%
  select(gene_id,
         gene_symbol,
         logFC_CRISPRoff = logFC, 
         padj_CRISPRoff = adj.P.Val) %>%
  inner_join(
    ko_results %>%
      select(gene_id,
             gene_symbol,
             logFC_KO = logFC,
             padj_KO = adj.P.Val),
    by = "gene_id"
  ) %>%
  # Use gene_symbol.x as the final gene_symbol (they should be the same)
  select(-gene_symbol.y) %>%
  rename(gene_symbol = gene_symbol.x) %>%
  mutate(
    significant = case_when(
      padj_CRISPRoff < 0.05 & padj_KO < 0.05 ~ "Both",
      padj_CRISPRoff < 0.05 ~ "CRISPRoff only",
      padj_KO < 0.05 ~ "KO only",
      TRUE ~ "Not significant"
    )
  )

# Create the comparison plot
comparison_plot <- ggplot(combined_results, 
                         aes(x = logFC_KO, 
                             y = logFC_CRISPRoff, 
                             color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = combined_results %>%
      mutate(
        max_abs_fc = pmax(abs(logFC_KO), abs(logFC_CRISPRoff))
      ) %>%
      filter(padj_KO < 0.05 | padj_CRISPRoff < 0.05) %>%
      slice_max(max_abs_fc, n = 30),
    aes(label = gene_symbol),
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(
      "Both" = "purple",
      "CRISPRoff only" = "blue",
      "KO only" = "red",
      "Not significant" = "grey"
    )
  ) +
  labs(
    title = "Comparison of MED12 perturbation methods",
    x = "log2 Fold Change (KO vs AAVS1)",
    y = "log2 Fold Change (CRISPRoff vs NTC)",
    color = "Differential Expression Significance"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Save the plot
ggsave("CRISPRoff_KO_MED12/plots_med12_dge/method_comparison.png", 
       comparison_plot, 
       width = 14,  # Wider to accommodate legend
       height = 10)  # Height for square plot area plus title

# Create summary statistics
summary_stats <- combined_results %>%
  summarise(
    total_genes = n(),
    sig_both = sum(significant == "Both"),
    sig_crispri = sum(significant == "CRISPRoff only"),
    sig_ko = sum(significant == "KO only"),
    not_sig = sum(significant == "Not significant")
  )

# Save the combined results
write_csv(combined_results, 
          "CRISPRoff_KO_MED12/results_med12_dge/combined_DE_results.csv")

# Print summary
print("Analysis complete. Results saved to:")
print("- CRISPRoff_KO_MED12/plots_med12_dge/method_comparison.png")
print("- CRISPRoff_KO_MED12/results_med12_dge/combined_DE_results.csv")
print("\nSummary statistics:")
print(summary_stats) 