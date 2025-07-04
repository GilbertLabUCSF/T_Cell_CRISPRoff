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

# Set directories and files
salmon_dir <- "results/star_salmon"
sample_sheet_file <- "samplesheet_nfcore.csv"
tx2gene_file <- "results/star_salmon/tx2gene.tsv"

# Create output directories
dir.create("results_med12_dge", showWarnings = FALSE)
dir.create("plots_med12_dge", showWarnings = FALSE)

# Read and process sample sheet
sample_sheet <- read_csv(sample_sheet_file) %>%
  mutate(
    Donor = str_extract(sample, "D[12]"),
    Condition = case_when(
      grepl("MED12", sample) ~ "MED12_KO",
      grepl("AAVS1", sample) ~ "AAVS1_KO"
    )
  ) %>%
  mutate(
    Donor = factor(Donor),
    Condition = factor(Condition, levels = c("AAVS1_KO", "MED12_KO"))
  )

# Get salmon directories and import data
all_dirs <- list.dirs(salmon_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- all_dirs[basename(all_dirs) %in% sample_sheet$sample]

files <- file.path(sample_dirs, "quant.sf")
names(files) <- basename(sample_dirs)

tx2gene <- read_tsv(tx2gene_file, col_names = c("tx", "gene"))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Ensure counts and sample sheet are in same order
sample_idx <- match(sample_sheet$sample, colnames(txi$counts))
txi$counts <- txi$counts[, sample_idx]

# Create DGEList object and filter low expression genes
dge <- DGEList(counts = txi$counts, samples = sample_sheet)

keep_expr <- filterByExpr(dge, group = dge$samples$Condition)
dge <- dge[keep_expr, keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# Export normalized counts
counts_df <- data.frame(
  gene_id = rownames(dge$counts),
  cpm(dge)
) %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db,
                             keys = gsub("\\..*$", "", gene_id),
                             keytype = "ENSEMBL", 
                             column = "SYMBOL",
                             multiVals = "first"))

write_csv(counts_df, file.path("results_med12_dge", "normalized_counts.csv"))

# Create normalized counts data frame for plotting
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

# Create expression plot for MED12
create_gene_expression_plot <- function(data, target_gene) {
  data_subset <- data %>%
    filter(gene_symbol == target_gene)
  
  ggplot(data_subset, aes(x = Condition, y = log2_cpm, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.6) +
    scale_fill_manual(values = c("AAVS1_KO" = "grey", "MED12_KO" = "darkred")) +
    labs(title = paste0(target_gene, " Expression"),
         y = "log2 CPM",
         x = "") +
    theme_pubr() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

med12_expression_plot <- create_gene_expression_plot(norm_counts_df, "MED12")
ggsave("plots_med12_dge/med12_expression.png", med12_expression_plot, 
       width = 5, height = 5)

# Differential expression analysis
design <- model.matrix(~ Donor + Condition, data = dge$samples)

v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# Get differential expression results
de_results <- topTable(fit, coef = "ConditionMED12_KO", number = Inf, sort.by = "P") %>%
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

write_csv(de_results, file.path("results_med12_dge", "DE_MED12_vs_AAVS1.csv"))

# Create volcano plots
create_volcano_plot <- function(df, highlight_gene = "MED12") {
  ggplot(df, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = ifelse(gene_symbol == highlight_gene, "Target",
                                 ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                      ifelse(logFC > 0, "Up", "Down"),
                                      "NotSig"))), 
               alpha = 0.6) +
    geom_text_repel(
      data = subset(df, gene_symbol == highlight_gene | 
                     (adj.P.Val < 0.05 & abs(logFC) > 2)),
      aes(label = gene_symbol),
      max.overlaps = 30
    ) +
    scale_color_manual(values = c("Up" = "red", 
                                 "Down" = "blue", 
                                 "NotSig" = "grey",
                                 "Target" = "orange"),
                      name = "Significance") +
    labs(title = "MED12 KO vs AAVS1 Control",
         x = "log2 Fold Change", 
         y = "-log10 P-value") +
    ggpubr::theme_pubr() +
    theme(legend.position = "bottom")
}

create_volcano_plot_adjusted <- function(df, highlight_gene = "MED12") {
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(gene_symbol == highlight_gene, "Target",
                                 ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                      ifelse(logFC > 0, "Up", "Down"),
                                      "NotSig"))), 
               alpha = 0.6) +
    geom_text_repel(
      data = subset(df, gene_symbol == highlight_gene | 
                     (adj.P.Val < 0.05 & abs(logFC) > 2)),
      aes(label = gene_symbol),
      max.overlaps = 30
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Up" = "red", 
                                 "Down" = "blue", 
                                 "NotSig" = "grey",
                                 "Target" = "orange"),
                      name = "Significance") +
    labs(title = "MED12 KO vs AAVS1 Control",
         x = "log2 Fold Change", 
         y = "-log10 adj.P-value") +
    ggpubr::theme_pubr() +
    theme(legend.position = "bottom")
}

volcano_plot <- create_volcano_plot(de_results)
volcano_plot_adjusted <- create_volcano_plot_adjusted(de_results)

ggsave("plots_med12_dge/volcano.png", volcano_plot, width = 8, height = 6)
ggsave("plots_med12_dge/volcano_adjusted.png", volcano_plot_adjusted, width = 8, height = 6)

# Print summary
print('Analysis complete. Results written to results_med12_dge/* and plots written to plots_med12_dge/*')