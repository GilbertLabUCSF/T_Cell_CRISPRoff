library(limma)
library(edgeR)
library(tximport)
library(org.Hs.eg.db)
library(tidyverse)
library(ComplexHeatmap)
library(ggrepel)
library(DT)
library(fgsea)
library(msigdbr)
library(patchwork)
library(ggpubr)

salmon_dir <- "results/star_salmon"
sample_sheet_file <- "samplesheet_CRISPRoff_DGE.csv"
tx2gene_file <- "results/star_salmon/tx2gene.tsv"

dir.create("results_dge", showWarnings = FALSE)
dir.create("plots_dge", showWarnings = FALSE)

sample_sheet <- read_csv(sample_sheet_file) %>%
  filter(Cas == "CRISPRoff") %>%  # Only keep CRISPRoff samples
  mutate(
    Donor = factor(Donor),
    Target = factor(Target, levels = c("NTC", "CD81", "CD55", "CD29", "CD151")),
    `Experiment Day` = factor(`Experiment Day`)
  ) %>%
  arrange(Target, Donor)

all_dirs <- list.dirs(salmon_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- all_dirs[basename(all_dirs) %in% sample_sheet$sample]

sample_sheet <- sample_sheet %>% 
  filter(sample %in% basename(sample_dirs))

files <- file.path(sample_dirs, "quant.sf")
names(files) <- basename(sample_dirs)

tx2gene <- read_tsv(tx2gene_file, col_names = c("tx", "gene"))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

sample_idx <- match(sample_sheet$sample, colnames(txi$counts))
txi$counts <- txi$counts[, sample_idx]

dge <- DGEList(counts = txi$counts, samples = sample_sheet)

keep_expr <- filterByExpr(dge, group = dge$samples$Target)
dge <- dge[keep_expr, keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

counts_df <- data.frame(
  gene_id = rownames(dge$counts),
  cpm(dge)
) %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db,
                             keys = gsub("\\..*$", "", gene_id),
                             keytype = "ENSEMBL", 
                             column = "SYMBOL",
                             multiVals = "first"))

write_csv(counts_df, file.path("results_dge", "normalized_counts.csv"))

norm_counts <- cpm(dge, log=TRUE)
norm_counts_df <- as.data.frame(norm_counts) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "log2_cpm") %>%
  left_join(sample_sheet, by = "sample") %>%
  mutate(gene_symbol = mapIds(
    org.Hs.eg.db,
    keys = gsub("\\..*$", "", gene_id),
    keytype = "ENSEMBL",
    column = "SYMBOL"
  ))

create_target_expression_plot <- function(data, target_gene) {
  # For CD29/ITGB1, need to handle both names
  data_subset <- data %>%
    filter(Target %in% c("NTC", ifelse(target_gene == "ITGB1", "CD29", target_gene))) %>%
    filter(gene_symbol == target_gene)
  
  # Print values being plotted
  print("Values being plotted:")
  print(data_subset %>% select(sample, Target, log2_cpm))
  
  target_label <- ifelse(target_gene == "ITGB1", "CD29", target_gene)
  
  ggplot(data_subset, aes(x = Target, y = log2_cpm, fill = Target)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.6) +
    scale_fill_manual(values = c("NTC" = "grey", target_label = "darkred")) +
    labs(title = paste0(target_gene, " Expression"),
         y = "log2 CPM", 
         x = "") +
    theme(legend.position = "none") +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ylim(-3, 8)
}

target_genes <- c("CD81", "CD55", "ITGB1", "CD151")

target_expression_plots <- lapply(target_genes, function(gene) {
  create_target_expression_plot(norm_counts_df, gene)
})

combined_target_plot <- wrap_plots(target_expression_plots, ncol = 2)
ggsave("plots_dge/target_expression_combined.png", combined_target_plot, 
       width = 10, height = 15)

# Design matrix: ~ Donor + Target
design <- model.matrix(~ Donor + Target, data = dge$samples)

v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

coefs <- colnames(design)
coefs_target <- coefs[grepl("^Target", coefs)]

all_results <- list()
for (coef_name in coefs_target) {
  target_gene <- gsub("Target", "", coef_name)
  tbl <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene_id") %>%
    mutate(
      gene_symbol = mapIds(
        org.Hs.eg.db,
        keys = gsub("\\..*$", "", gene_id),
        keytype = "ENSEMBL",
        column = "SYMBOL"
      )
    ) %>% 
    select(gene_id, gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
  
  out_file <- file.path("results_dge", paste0("DE_", target_gene, "_vs_NTC.csv"))
  write_csv(tbl, out_file)
  
  all_results[[target_gene]] <- tbl
}

create_volcano_plot <- function(df, target_gene) {
  ggplot(df, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = ifelse(is.na(gene_symbol), 
                                  "NA", 
                                  ifelse(gene_symbol == target_gene, "Target",
                                         ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                                ifelse(logFC > 0, "Up", "Down"),
                                                ifelse(is.na(adj.P.Val), "NA", "NotSig"))))), 
               alpha = 0.6) +
    geom_text_repel(
      data = subset(df, gene_symbol == target_gene),
      aes(label = gene_symbol),
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c("Up" = "red", 
                                 "Down" = "blue", 
                                 "NotSig" = "grey30",
                                 "Target" = "orange",
                                 "NA" = "grey30"),
                      name = "Significance") +
    labs(title = paste0(target_gene, " vs NTC"),
         x = "log2 Fold Change", 
         y = "-log10 P-value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(t = 10, b = 10))
}

create_volcano_plot_adjusted <- function(df, target_gene) {
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(gene_symbol == target_gene, "Target",
                                 ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                      ifelse(logFC > 0, "Up", "Down"),
                                      "NotSig"))), 
               alpha = 0.6) +
    geom_text_repel(
      data = subset(df, gene_symbol == target_gene),
      aes(label = gene_symbol),
      max.overlaps = Inf
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Up" = "red", 
                                 "Down" = "blue", 
                                 "NotSig" = "grey30",
                                 "Target" = "orange",
                                 "NA" = "grey30"),
                      name = "Significance") +
    labs(title = paste0(target_gene, " vs NTC"),
         x = "log2 Fold Change", 
         y = "-log10 adj.P-value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(t = 10, b = 10))
}

volcano_plots <- list(
  create_volcano_plot(all_results[["CD81"]], "CD81"),
  create_volcano_plot(all_results[["CD55"]], "CD55"),
  create_volcano_plot(all_results[["CD29"]], "CD29"),
  create_volcano_plot(all_results[["CD151"]], "CD151")
)

combined_volcano_plot <- wrap_plots(volcano_plots, ncol = 2) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

volcano_plots_adjusted <- list(
  create_volcano_plot_adjusted(all_results[["CD81"]], "CD81"),
  create_volcano_plot_adjusted(all_results[["CD55"]], "CD55"),
  create_volcano_plot_adjusted(all_results[["CD29"]], "CD29"),
  create_volcano_plot_adjusted(all_results[["CD151"]], "CD151")
)

combined_volcano_plot_adjusted <- wrap_plots(volcano_plots_adjusted, ncol = 2) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


ggsave("plots_dge/volcano_combined.png", combined_volcano_plot, width = 12, height = 6)
ggsave("plots_dge/volcano_combined_adjusted.png", combined_volcano_plot_adjusted, width = 12, height = 6)
print('Analysis complete. Results written to results_dge/* and plots written to plots_dge/*') 