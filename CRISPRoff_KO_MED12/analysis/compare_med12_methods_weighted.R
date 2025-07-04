library(tidyverse)
library(ggrepel)
library(ggpubr)

# Read the DE results
crispri_results <- read_csv("CRISPRoff/results_dge/DE_targetMED12_vs_NTC.csv")
ko_results <- read_csv("CRISPRoff_KO_MED12/results_med12_dge/DE_MED12_vs_AAVS1.csv")

# Combine the results with weighted scores
combined_results <- crispri_results %>%
  select(gene_id,
         gene_symbol,
         logFC_CRISPRoff = logFC, 
         padj_CRISPRoff = adj.P.Val,
         pvalue_CRISPRoff = P.Value) %>%
  inner_join(
    ko_results %>%
      select(gene_id,
             gene_symbol,
             logFC_KO = logFC,
             padj_KO = adj.P.Val,
             pvalue_KO = P.Value),
    by = "gene_id"
  ) %>%
  select(-gene_symbol.y) %>%
  rename(gene_symbol = gene_symbol.x) %>%
  mutate(
    # Calculate weighted scores (LFC * -log10(p-value))
    weighted_score_CRISPRoff = logFC_CRISPRoff * -log10(pvalue_CRISPRoff),
    weighted_score_KO = logFC_KO * -log10(pvalue_KO),
    significant = case_when(
      padj_CRISPRoff < 0.05 & padj_KO < 0.05 ~ "Both",
      padj_CRISPRoff < 0.05 ~ "CRISPRoff only",
      padj_KO < 0.05 ~ "KO only",
      TRUE ~ "Not significant"
    )
  )

# Create the comparison plot with weighted scores
comparison_plot_weighted <- ggplot(combined_results, 
                                 aes(x = weighted_score_KO, 
                                     y = weighted_score_CRISPRoff, 
                                     color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = combined_results %>%
      mutate(
        max_abs_weighted = pmax(abs(weighted_score_KO), 
                               abs(weighted_score_CRISPRoff))
      ) %>%
      filter(padj_KO < 0.05 | padj_CRISPRoff < 0.05) %>%
      slice_max(max_abs_weighted, n = 30),
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
    subtitle = "Weighted by effect size and significance",
    x = "Weighted Score (KO vs AAVS1)\nlog2FC * -log10(p-value)",
    y = "Weighted Score (CRISPRoff vs NTC)\nlog2FC * -log10(p-value)",
    color = "Differential Expression Significance"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the weighted plot
ggsave("CRISPRoff_KO_MED12/plots_med12_dge/method_comparison_weighted.png", 
       comparison_plot_weighted, 
       width = 14,  # Wider to accommodate legend
       height = 10)  # Height for square plot area plus title

# Save the combined results with weighted scores
write_csv(combined_results, 
          "CRISPRoff_KO_MED12/results_med12_dge/combined_DE_results_weighted.csv")

# Print summary
print("Analysis complete. Results saved to:")
print("- CRISPRoff_KO_MED12/plots_med12_dge/method_comparison_weighted.png")
print("- CRISPRoff_KO_MED12/results_med12_dge/combined_DE_results_weighted.csv")

# Print summary of weighted scores
summary_weighted <- combined_results %>%
  summarise(
    total_genes = n(),
    mean_weight_CRISPRoff = mean(abs(weighted_score_CRISPRoff), na.rm = TRUE),
    mean_weight_KO = mean(abs(weighted_score_KO), na.rm = TRUE),
    max_weight_CRISPRoff = max(abs(weighted_score_CRISPRoff), na.rm = TRUE),
    max_weight_KO = max(abs(weighted_score_KO), na.rm = TRUE)
  )

print("\nSummary of weighted scores:")
print(summary_weighted) 