# Load Libraries
library(limma)
library(edgeR)
library(tximport)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(ComplexHeatmap)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(patchwork)
library(ggpubr)
library(magrittr)
library(ggVennDiagram)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(circlize)
library(MASS)

# Configuration
EXPT_NAME                <- "Co065"
NTC_GENE_SYMBOL          <- "AAVS1"
DONOR_TO_EXCLUDE         <- ""           # leave empty for none
P_VALUE_CUTOFF           <- 0.05
LOGFC_CUTOFF             <- 1
MIN_GENES_FOR_ENRICHMENT <- 10

# targets studied in this experiment
targets <- c("FAS", "PTPN2", "RC3H1", "SUV39H1", "RASA2", "MED12")

# Read Input Data
sample_sheet    <- read_csv("rna_seq_meta_key.csv")
txi             <- readRDS(file.path("data", "tximport_object.rds"))
dge             <- readRDS(file.path("data", "dge_object.rds"))
norm_counts_csv <- read.csv(file.path("data", "normalized_counts.csv"))

# Initial DGE Object Filtering
if (!"sample" %in% colnames(dge$samples)) {
  if (all(rownames(dge$samples) %in% colnames(dge$counts))) {
    warning("dge$samples lacks a 'sample' column; using rownames.")
    dge$samples$sample <- rownames(dge$samples)
  } else stop("dge$samples must contain a 'sample' column.")
}

required_cols <- c("Experiment.ID", "Target", "Cas", "Donor")
missing_cols  <- setdiff(required_cols, colnames(dge$samples))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

cat("Original samples:", nrow(dge$samples), "\n")

dge$samples <- dge$samples %>% filter(Experiment.ID == EXPT_NAME)
if (DONOR_TO_EXCLUDE != "")
  dge$samples <- dge$samples %>% filter(Donor != DONOR_TO_EXCLUDE)

if (nrow(dge$samples) == 0) stop("No samples left after filtering.\n")

dge$counts         <- dge$counts[, dge$samples$sample]
dge$samples$Target <- toupper(dge$samples$Target)
stopifnot(all(colnames(dge$counts) == dge$samples$sample))

# Expression Boxplot Analysis

# Gene Filtering & Design
keep.exprs <- filterByExpr(
  dge,
  group            = interaction(dge$samples$Target, dge$samples$Cas),
  min.count        = 10,
  min.total.count  = 15
)
dge <- dge[keep.exprs, , keep.lib.sizes = FALSE]
if (nrow(dge) == 0) stop("All genes filtered out by filterByExpr.\n")

dge$samples$condition <- factor(paste0(dge$samples$Cas, "_", dge$samples$Target))
dge$samples$Donor     <- factor(gsub(" ", "_", dge$samples$Donor))

design <- model.matrix(~ 0 + condition + Donor, data = dge$samples)
colnames(design) <- make.names(colnames(design))

# Voom + Limma Analysis
v <- voom(dge, design, plot = FALSE)

contr.matrix <- makeContrasts(
  Fas_KO     = Cas9_FAS        - Cas9_AAVS1,
  Fas_KD     = CRISPRoff_FAS   - CRISPRoff_AAVS1,
  PTPN2_KO   = Cas9_PTPN2      - Cas9_AAVS1,
  PTPN2_KD   = CRISPRoff_PTPN2 - CRISPRoff_AAVS1,
  RC3H1_KO   = Cas9_RC3H1      - Cas9_AAVS1,
  RC3H1_KD   = CRISPRoff_RC3H1 - CRISPRoff_AAVS1,
  SUV39H1_KO = Cas9_SUV39H1    - Cas9_AAVS1,
  SUV39H1_KD = CRISPRoff_SUV39H1-CRISPRoff_AAVS1,
  RASA2_KO   = Cas9_RASA2      - Cas9_AAVS1,
  RASA2_KD   = CRISPRoff_RASA2 - CRISPRoff_AAVS1,
  MED12_KO   = Cas9_MED12      - Cas9_AAVS1,
  MED12_KD   = CRISPRoff_MED12 - CRISPRoff_AAVS1,
  # double-differential
  Fas_KDvsKO      = (CRISPRoff_FAS     - CRISPRoff_AAVS1) - (Cas9_FAS     - Cas9_AAVS1),
  PTPN2_KDvsKO    = (CRISPRoff_PTPN2   - CRISPRoff_AAVS1) - (Cas9_PTPN2   - Cas9_AAVS1),
  RC3H1_KDvsKO    = (CRISPRoff_RC3H1   - CRISPRoff_AAVS1) - (Cas9_RC3H1   - Cas9_AAVS1),
  SUV39H1_KDvsKO  = (CRISPRoff_SUV39H1 - CRISPRoff_AAVS1) - (Cas9_SUV39H1 - Cas9_AAVS1),
  RASA2_KDvsKO    = (CRISPRoff_RASA2   - CRISPRoff_AAVS1) - (Cas9_RASA2   - Cas9_AAVS1),
  MED12_KDvsKO    = (CRISPRoff_MED12   - CRISPRoff_AAVS1) - (Cas9_MED12   - Cas9_AAVS1),
  levels = design
)

fit <- lmFit(v, design) |>
       contrasts.fit(contr.matrix) |>
       eBayes(robust = TRUE)

# Save Top-Tables
dir.create(file.path("results", EXPT_NAME), showWarnings = FALSE, recursive = TRUE)
all_results_topTable <- list()

for (coef_name in colnames(fit$contrasts)) {
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") |>
    rownames_to_column("gene_id") |>
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = gsub("\\..*$", "", gene_id),
                                keytype = "ENSEMBL", column = "SYMBOL",
                                multiVals = "first"),
           entrez      = mapIds(org.Hs.eg.db,
                                keys = gsub("\\..*$", "", gene_id),
                                keytype = "ENSEMBL", column = "ENTREZID",
                                multiVals = "first"))
  all_results_topTable[[coef_name]] <- tt
  write_csv(tt, file.path("results", EXPT_NAME,
                          paste0("DE_topTable_", EXPT_NAME, "_", coef_name, ".csv")))
}

# Pre-compute KO∩KD significant genes
get_sig <- function(df) {
  if (is.null(df)) return(character(0))
  df |>
    filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) |>
    pull(gene_symbol) |>
    na.omit() |>
    unique()
}

common_sig_by_target <- setNames(vector("list", length(targets)), targets)
for (tg in targets) {
  common_sig_by_target[[tg]] <- intersect(
    get_sig(all_results_topTable[[paste0(tg, "_KO")]]),
    get_sig(all_results_topTable[[paste0(tg, "_KD")]])
  )
}

# Enhanced Volcano Plot
create_volcano_plot <- function(df_results, contrast_id,
                                top_effect_n = 5,
                                top_p_n      = 5,
                                custom_genes = NULL,
                                common_genes = NULL,
                                logfc_cutoff = LOGFC_CUTOFF,
                                p_cutoff     = P_VALUE_CUTOFF) {

  df_results <- df_results |>
    mutate(Significance = case_when(
      adj.P.Val < p_cutoff & logFC >  logfc_cutoff ~ "Up",
      adj.P.Val < p_cutoff & logFC < -logfc_cutoff ~ "Down",
      TRUE                                          ~ "NotSig"
    ))

  # Limit candidate pool to genes present in common_genes (if supplied)
  pool_df <- if (!is.null(common_genes) && length(common_genes) > 0) {
               df_results %>% filter(gene_symbol %in% common_genes)
             } else df_results

  top_eff <- pool_df |> arrange(-abs(logFC)) |> slice_head(n = top_effect_n) |> pull(gene_symbol)
  top_pv  <- pool_df |> arrange(P.Value)      |> slice_head(n = top_p_n)     |> pull(gene_symbol)

  to_label <- union(union(top_eff, top_pv), custom_genes) |> na.omit()

  ggplot(df_results, aes(logFC, -log10(P.Value))) +
    geom_hline(yintercept = -log10(0.01), linetype = "dotted", colour = "grey70") +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
               linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = Significance), alpha = 0.75, size = 1.5) +
    scale_colour_manual(values = c(Up = "red", Down = "blue", NotSig = "grey60")) +
    ggrepel::geom_text_repel(
      data = df_results %>% filter(gene_symbol %in% to_label),
      aes(label = gene_symbol),
      box.padding = 0.4, point.padding = 0.15,
      max.overlaps = Inf, size = 3
    ) +
    labs(title = contrast_id,
         x = "log2 Fold Change",
         y = "-log10 P-value") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
}

# Build & Save All Volcano Plots
volcano_plot_params <- unlist(lapply(targets, function(tg) {
  list(
    list(contrast = paste0(tg, "_KO"),     tg = tg),
    list(contrast = paste0(tg, "_KD"),     tg = tg),
    list(contrast = paste0(tg, "_KDvsKO"), tg = tg)  # include DD if desired
  )
}), recursive = FALSE)

volcano_plots <- lapply(volcano_plot_params, function(p) {
  df <- all_results_topTable[[p$contrast]]
  if (is.null(df)) return(NULL)

  create_volcano_plot(
    df_results  = df,
    contrast_id = p$contrast,
    top_effect_n = 5,
    top_p_n      = 5,
    custom_genes = NULL,
    common_genes = common_sig_by_target[[p$tg]]
  )
})

volcano_plots <- Filter(Negate(is.null), volcano_plots)

if (length(volcano_plots) > 0) {
  dir.create(file.path("plots", EXPT_NAME), showWarnings = FALSE, recursive = TRUE)

  combined <- wrap_plots(volcano_plots, ncol = 3, guides = "collect") &
              theme(legend.position = "bottom")

  ggsave(file.path("plots", EXPT_NAME,
                   paste0(EXPT_NAME, "_combined_volcano.png")),
         combined,
         width = 15, height = 5 * ceiling(length(volcano_plots) / 3))
}

# Gene Ontology and Pathway Enrichment Analysis

# Enhanced ORA Functions
perform_go_enrichment <- function(gene_list, universe_genes, ont = "BP", pval_cutoff = 0.05) {
  ego <- enrichGO(gene = gene_list,
                  universe = universe_genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = pval_cutoff,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  return(ego)
}

perform_kegg_enrichment <- function(gene_list, universe_genes, pval_cutoff = 0.05) {
  ekegg <- enrichKEGG(gene = gene_list,
                      universe = universe_genes,
                      organism = 'hsa',
                      pAdjustMethod = "BH",
                      pvalueCutoff = pval_cutoff,
                      qvalueCutoff = 0.2)
  return(ekegg)
}

# Get universe of all tested genes
universe_genes <- unique(na.omit(unlist(lapply(all_results_topTable, function(x) x$entrez))))

# Create enrichment results directories
if (!dir.exists("results_Co065/go")) dir.create("results_Co065/go", recursive = TRUE)
if (!dir.exists("results_Co065/enrichment_comparison")) dir.create("results_Co065/enrichment_comparison", recursive = TRUE)

# Perform enrichment for each contrast
go_results <- list()
kegg_results <- list()
enrichment_summary <- list()

for (contrast_name in names(all_results_topTable)) {
  cat("\nPerforming enrichment analysis for:", contrast_name, "\n")
  
  # Get significant genes
  sig_genes <- all_results_topTable[[contrast_name]] %>%
    filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) %>%
    pull(entrez) %>%
    na.omit() %>%
    unique()
  
  if (length(sig_genes) >= MIN_GENES_FOR_ENRICHMENT) {
    # GO enrichment
    ego <- perform_go_enrichment(sig_genes, universe_genes)
    go_results[[contrast_name]] <- ego
    
    # KEGG enrichment
    ekegg <- perform_kegg_enrichment(sig_genes, universe_genes)
    kegg_results[[contrast_name]] <- ekegg
    
    # Store summary for comparison
    if (nrow(as.data.frame(ego)) > 0) {
      enrichment_summary[[contrast_name]] <- as.data.frame(ego) %>%
        dplyr::select(ID, Description, p.adjust, Count) %>%
        mutate(Contrast = contrast_name)
    }
    
    # Save results
    if (nrow(as.data.frame(ego)) > 0) {
      write_csv(as.data.frame(ego), file.path("results_Co065/go", paste0("GO_", EXPT_NAME, "_", contrast_name, ".csv")))
      
      # Create individual plots
      p1 <- dotplot(ego, showCategory = 20) + 
        ggtitle(paste("GO Enrichment:", contrast_name)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave(file.path("plots_dge", paste0(EXPT_NAME, "_GO_", contrast_name, ".png")), 
             p1, width = 8, height = 10)
    }
    
    if (nrow(as.data.frame(ekegg)) > 0) {
      if (!dir.exists("results_go")) dir.create("results_go", recursive = TRUE)
      write_csv(as.data.frame(ekegg), file.path("results_go", paste0("KEGG_", EXPT_NAME, "_", contrast_name, ".csv")))
    }
  } else {
    cat("Too few significant genes for enrichment in", contrast_name, "\n")
  }
}

# Comprehensive Comparison of Enriched Terms
# Create comparison matrix for GO terms
create_enrichment_comparison_matrix <- function(go_results_list, top_n = 30) {
  # Collect all unique GO terms
  all_go_terms <- c()
  go_data_list <- list()
  
  for (contrast_name in names(go_results_list)) {
    if (!is.null(go_results_list[[contrast_name]])) {
      df <- as.data.frame(go_results_list[[contrast_name]])
      if (nrow(df) > 0) {
        go_data_list[[contrast_name]] <- df
        all_go_terms <- c(all_go_terms, df$Description)
      }
    }
  }
  
  all_go_terms <- unique(all_go_terms)
  
  if (length(all_go_terms) == 0) return(NULL)
  
  # Create matrix of -log10(p.adjust)
  mat <- matrix(0, nrow = length(all_go_terms), ncol = length(go_data_list))
  rownames(mat) <- all_go_terms
  colnames(mat) <- names(go_data_list)
  
  for (contrast in names(go_data_list)) {
    df <- go_data_list[[contrast]]
    idx <- match(df$Description, all_go_terms)
    mat[idx, contrast] <- -log10(df$p.adjust)
  }
  
  # Select top terms by max significance across contrasts
  if (nrow(mat) > top_n) {
    top_idx <- order(apply(mat, 1, max), decreasing = TRUE)[1:top_n]
    mat <- mat[top_idx, , drop = FALSE]
  }
  
  return(mat)
}

# Create comparison heatmap
comparison_mat <- create_enrichment_comparison_matrix(go_results, top_n = 40)

if (!is.null(comparison_mat)) {
  # Order columns by similarity
  col_order <- hclust(dist(t(comparison_mat)))$order
  comparison_mat <- comparison_mat[, col_order]
  
  # Create heatmap
  col_fun <- colorRamp2(c(0, 1.3, 3, 5), c("white", "lightblue", "blue", "darkblue"))
  
  ht <- Heatmap(comparison_mat,
                name = "-log10(adj.P)",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                row_names_max_width = unit(10, "cm"),
                column_names_rot = 45,
                column_title = "GO Term Enrichment Across All Contrasts",
                column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 10),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if (comparison_mat[i, j] > 1.3) {
                    grid.text("*", x, y, gp = gpar(fontsize = 8))
                  }
                })
  
  if (!dir.exists("results_enrichment_comparison")) dir.create("results_enrichment_comparison", recursive = TRUE)
  png(file.path("results_enrichment_comparison", paste0(EXPT_NAME, "_GO_comparison_heatmap.png")), 
      width = 12, height = 14, units = "in", res = 300)
  draw(ht)
  dev.off()
}

# Find Common Enriched Terms
find_common_enriched_terms <- function(go_results_list, p_cutoff = 0.05) {
  term_counts <- list()
  term_details <- list()
  
  for (contrast in names(go_results_list)) {
    if (!is.null(go_results_list[[contrast]])) {
      df <- as.data.frame(go_results_list[[contrast]])
      sig_terms <- df %>% filter(p.adjust < p_cutoff)
      
      for (i in seq_len(nrow(sig_terms))) {
        term <- sig_terms$Description[i]
        if (!(term %in% names(term_counts))) {
          term_counts[[term]] <- 0
          term_details[[term]] <- list(contrasts = c(), p_values = c())
        }
        term_counts[[term]] <- term_counts[[term]] + 1
        term_details[[term]]$contrasts <- c(term_details[[term]]$contrasts, contrast)
        term_details[[term]]$p_values <- c(term_details[[term]]$p_values, sig_terms$p.adjust[i])
      }
    }
  }
  
  # Create summary dataframe
  common_terms_df <- data.frame(
    Term = names(term_counts),
    Count = unlist(term_counts),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Count)) %>%
    filter(Count > 1)
  
  # Add details
  common_terms_df$Contrasts <- sapply(common_terms_df$Term, function(x) 
    paste(term_details[[x]]$contrasts, collapse = ", "))
  common_terms_df$Mean_adj_P <- sapply(common_terms_df$Term, function(x) 
    mean(term_details[[x]]$p_values))
  
  return(common_terms_df)
}

# Find common terms
common_go_terms <- find_common_enriched_terms(go_results)
if (nrow(common_go_terms) > 0) {
  write_csv(common_go_terms, file.path("results_enrichment_comparison", 
                                       paste0(EXPT_NAME, "_common_GO_terms.csv")))
  
  # Visualize common terms
  p_common <- ggplot(common_go_terms %>% head(20), 
                     aes(x = reorder(Term, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "GO Terms Enriched in Multiple Contrasts",
         x = "GO Term",
         y = "Number of Contrasts") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))
  
  ggsave(file.path("results_enrichment_comparison", 
                   paste0(EXPT_NAME, "_common_GO_terms_barplot.png")), 
         p_common, width = 10, height = 8)
}

# Comparison between specific contrast groups
# KO vs KD comparison for each target
targets <- c("FAS", "PTPN2", "RC3H1", "SUV39H1", "RASA2", "MED12")

for (target in targets) {
  ko_contrast <- paste0(target, "_KO")
  kd_contrast <- paste0(target, "_KD")
  
  if (ko_contrast %in% names(go_results) && kd_contrast %in% names(go_results)) {
    ko_df <- as.data.frame(go_results[[ko_contrast]])
    kd_df <- as.data.frame(go_results[[kd_contrast]])
    
    if (nrow(ko_df) > 0 && nrow(kd_df) > 0) {
      # Find overlapping terms
      common_terms <- intersect(ko_df$ID, kd_df$ID)
      ko_only <- setdiff(ko_df$ID, kd_df$ID)
      kd_only <- setdiff(kd_df$ID, ko_df$ID)
      
      # Create detailed comparison
      comparison_df <- data.frame(
        Category = c(rep("Common", length(common_terms)),
                     rep("KO_only", length(ko_only)),
                     rep("KD_only", length(kd_only))),
        GO_ID = c(common_terms, ko_only, kd_only),
        stringsAsFactors = FALSE
      )
      
      # Add descriptions and p-values
      comparison_df <- comparison_df %>%
        left_join(bind_rows(
          ko_df %>% dplyr::select(ID, Description, p.adjust) %>% mutate(Source = "KO"),
          kd_df %>% dplyr::select(ID, Description, p.adjust) %>% mutate(Source = "KD")
        ) %>% 
          group_by(ID, Description) %>% 
          summarise(min_p.adjust = min(p.adjust), .groups = 'drop'),
        by = c("GO_ID" = "ID"))
      
      write_csv(comparison_df, file.path("results_enrichment_comparison", 
                                         paste0(EXPT_NAME, "_", target, "_KOvsKD_GO_comparison.csv")))
      
      # Create Venn diagram
      venn_data <- list(
        KO = ko_df$Description[ko_df$p.adjust < 0.05],
        KD = kd_df$Description[kd_df$p.adjust < 0.05]
      )
      
      venn_plot <- ggVennDiagram(venn_data, label_alpha = 0) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = paste("GO Terms:", target, "KO vs KD")) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave(file.path("results_enrichment_comparison", 
                       paste0(EXPT_NAME, "_", target, "_KOvsKD_GO_venn.png")), 
             venn_plot, width = 6, height = 5)
    }
  }
}

# MSigDB Gene Set Analysis
# Use MSigDB hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ncbi_gene)

hallmark_results <- list()
for (contrast_name in names(all_results_topTable)) {
  sig_genes <- all_results_topTable[[contrast_name]] %>%
    filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) %>%
    pull(entrez) %>%
    na.omit() %>%
    unique()
  
  if (length(sig_genes) >= MIN_GENES_FOR_ENRICHMENT) {
    hallmark_enrich <- enricher(gene = as.character(sig_genes),
                                universe = as.character(universe_genes),
                                TERM2GENE = hallmark_sets,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
    
    if (!is.null(hallmark_enrich) && nrow(as.data.frame(hallmark_enrich)) > 0) {
      hallmark_results[[contrast_name]] <- hallmark_enrich
      if (!dir.exists("results_go")) dir.create("results_go", recursive = TRUE)
      write_csv(as.data.frame(hallmark_enrich), 
                file.path("results_go", paste0("Hallmark_", EXPT_NAME, "_", contrast_name, ".csv")))
    }
  }
}

# Create Overall Enrichment Summary Report
create_enrichment_summary_report <- function() {
  report_lines <- c(
    paste("# Enrichment Analysis Summary for", EXPT_NAME),
    paste("Generated on:", Sys.Date()),
    "",
    "## Overview",
    paste("Total contrasts analyzed:", length(all_results_topTable)),
    "",
    "## Common GO Terms Across Multiple Contrasts"
  )
  
  if (exists("common_go_terms") && nrow(common_go_terms) > 0) {
    report_lines <- c(report_lines,
                      paste("Found", nrow(common_go_terms), "GO terms enriched in multiple contrasts"),
                      "",
                      "### Top 10 Most Common GO Terms:",
                      sapply(1:min(10, nrow(common_go_terms)), function(i) {
                        paste0(i, ". ", common_go_terms$Term[i], " (", common_go_terms$Count[i], " contrasts)")
                      }),
                      "")
  }
  
  # Add summary for each target
  report_lines <- c(report_lines, "", "## Target-Specific Summary")
  
  for (target in targets) {
    ko_name <- paste0(target, "_KO")
    kd_name <- paste0(target, "_KD")
    dd_name <- paste0(target, "_KDvsKO")
    
    report_lines <- c(report_lines, "", paste("###", target))
    
    # Add significant gene counts
    if (ko_name %in% names(all_results_topTable)) {
      ko_sig <- sum(all_results_topTable[[ko_name]]$adj.P.Val < P_VALUE_CUTOFF & 
                      abs(all_results_topTable[[ko_name]]$logFC) > LOGFC_CUTOFF, na.rm = TRUE)
      report_lines <- c(report_lines, paste("- KO significant genes:", ko_sig))
    }
    
    if (kd_name %in% names(all_results_topTable)) {
      kd_sig <- sum(all_results_topTable[[kd_name]]$adj.P.Val < P_VALUE_CUTOFF & 
                      abs(all_results_topTable[[kd_name]]$logFC) > LOGFC_CUTOFF, na.rm = TRUE)
      report_lines <- c(report_lines, paste("- KD significant genes:", kd_sig))
    }
    
    if (dd_name %in% names(all_results_topTable)) {
      dd_sig <- sum(all_results_topTable[[dd_name]]$adj.P.Val < P_VALUE_CUTOFF & 
                      abs(all_results_topTable[[dd_name]]$logFC) > LOGFC_CUTOFF, na.rm = TRUE)
      report_lines <- c(report_lines, paste("- KD vs KO differential genes:", dd_sig))
    }
  }
  
  writeLines(report_lines, file.path("results_enrichment_comparison", 
                                     paste0(EXPT_NAME, "_enrichment_summary_report.txt")))
}

create_enrichment_summary_report()

# Interactive HTML Report for Enrichment Results
create_interactive_enrichment_table <- function() {
  if (length(enrichment_summary) > 0) {
    all_enrichments <- bind_rows(enrichment_summary)
    
    # Create a matrix showing which contrasts have which terms
    term_contrast_matrix <- all_enrichments %>%
      mutate(Present = 1) %>%
      pivot_wider(names_from = Contrast, values_from = Present, values_fill = 0, 
                  id_cols = c(ID, Description))
    
    write_csv(term_contrast_matrix, 
              file.path("results_enrichment_comparison", 
                        paste0(EXPT_NAME, "_enrichment_presence_matrix.csv")))
  }
}

create_interactive_enrichment_table()

# Volcano Plots (including double differential)
create_volcano_plot <- function(df_results, contrast_id, gene_to_highlight = NULL) {
  df_results <- df_results %>%
    mutate(Significance = case_when(
      !is.na(gene_symbol) & gene_symbol == gene_to_highlight ~ "TargetGene",
      adj.P.Val < P_VALUE_CUTOFF & abs(logFC) > LOGFC_CUTOFF & logFC > 0 ~ "Up",
      adj.P.Val < P_VALUE_CUTOFF & abs(logFC) > LOGFC_CUTOFF & logFC < 0 ~ "Down",
      TRUE ~ "NotSig"
    )) %>%
    mutate(Significance = factor(Significance,
                                 levels = c("Up","Down","NotSig","TargetGene")))
  
  p <- ggplot(df_results, aes(logFC, -log10(P.Value), colour = Significance)) +
    geom_hline(yintercept = -log10(0.01), linetype = "dotted", colour = "grey70") +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_colour_manual(
      breaks = levels(df_results$Significance),
      values = c(Up = "red", Down = "blue", NotSig = "grey60", TargetGene = "orange"),
      drop = FALSE,
      name = "Significance"
    ) +
    labs(title = contrast_id,
         x = "log2 Fold Change",
         y = "-log10 P-value") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(gene_to_highlight) && gene_to_highlight %in% df_results$gene_symbol) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(df_results, Significance == "TargetGene"),
      aes(label = gene_symbol),
      max.overlaps = Inf, box.padding = 0.5, point.padding = 0.2,
      colour = "black", fontface = "bold"
    )
  }
  
  return(p)
}

# Create volcano plots for all contrasts including double differential
volcano_plot_params <- list(
  list(contrast = "Fas_KO", highlight_gene = "FAS"),
  list(contrast = "Fas_KD", highlight_gene = "FAS"),
  list(contrast = "Fas_KDvsKO", highlight_gene = "FAS"),
  list(contrast = "PTPN2_KO", highlight_gene = "PTPN2"),
  list(contrast = "PTPN2_KD", highlight_gene = "PTPN2"),
  list(contrast = "PTPN2_KDvsKO", highlight_gene = "PTPN2"),
  list(contrast = "RC3H1_KO", highlight_gene = "RC3H1"),
  list(contrast = "RC3H1_KD", highlight_gene = "RC3H1"),
  list(contrast = "RC3H1_KDvsKO", highlight_gene = "RC3H1"),
  list(contrast = "SUV39H1_KO", highlight_gene = "SUV39H1"),
  list(contrast = "SUV39H1_KD", highlight_gene = "SUV39H1"),
  list(contrast = "SUV39H1_KDvsKO", highlight_gene = "SUV39H1"),
  list(contrast = "RASA2_KO", highlight_gene = "RASA2"),
  list(contrast = "RASA2_KD", highlight_gene = "RASA2"),
  list(contrast = "RASA2_KDvsKO", highlight_gene = "RASA2"),
  list(contrast = "MED12_KO", highlight_gene = "MED12"),
  list(contrast = "MED12_KD", highlight_gene = "MED12"),
  list(contrast = "MED12_KDvsKO", highlight_gene = "MED12")
)

volcano_plots_list <- lapply(volcano_plot_params, function(p) {
  if (p$contrast %in% names(all_results_topTable)) {
    create_volcano_plot(all_results_topTable[[p$contrast]], p$contrast, p$highlight_gene)
  } else { NULL }
})
volcano_plots_list <- Filter(Negate(is.null), volcano_plots_list)

if (length(volcano_plots_list) > 0) {
  combined_volcano_plot <- wrap_plots(volcano_plots_list, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path("plots_dge", paste0(EXPT_NAME, "_combined_volcano_all_contrasts.png")), 
         combined_volcano_plot,
         width = 15, height = 5 * ceiling(length(volcano_plots_list)/3))
}

# Correlation Plots
create_corr_plot_robust_rlm <- function(df1, df2, gene_highlight, labA, labB,
                                        colA = "red", colB = "blue",
                                        colBoth = "purple", colNS = "grey70",
                                        rlm_method_psi = "psi.bisquare") {
  
  if (!"gene_symbol" %in% names(df1)) df1 <- df1 %>% mutate(gene_symbol = gene_id)
  if (!"gene_symbol" %in% names(df2)) df2 <- df2 %>% mutate(gene_symbol = gene_id)
  
  df1_sub <- df1 %>% dplyr::select(gene_id, gene_symbol, logFC_df1 = logFC, adj.P.Val_df1 = adj.P.Val)
  df2_sub <- df2 %>% dplyr::select(gene_id, gene_symbol, logFC_df2 = logFC, adj.P.Val_df2 = adj.P.Val)
  
  dge_tab <- inner_join(df1_sub, df2_sub, by = "gene_id") %>%
    mutate(gene_symbol = ifelse(!is.na(gene_symbol.x), gene_symbol.x, gene_symbol.y)) %>%
    dplyr::select(-any_of(c("gene_symbol.x", "gene_symbol.y"))) %>%
    rename(logFC_A = logFC_df1, adj.P.Val_A = adj.P.Val_df1, logFC_B = logFC_df2, adj.P.Val_B = adj.P.Val_df2) %>%
    mutate(
      sigA = adj.P.Val_A < P_VALUE_CUTOFF,
      sigB = adj.P.Val_B < P_VALUE_CUTOFF,
      significant_cat = case_when(
        sigA & sigB ~ "Both",
        sigA & !sigB ~ labA,
        !sigA & sigB ~ labB,
        TRUE ~ "Not significant"
      ),
      significant_cat = factor(significant_cat, levels = c("Both", labA, labB, "Not significant"))
    )
  
  robust_line_data <- NULL
  if (nrow(dge_tab) >= 5) {
    tryCatch({
      rlm_fit <- MASS::rlm(logFC_B ~ logFC_A, data = dge_tab, psi = rlm_method_psi, method="M", maxit=100)
      x_range <- range(dge_tab$logFC_A, na.rm = TRUE)
      x_seq <- seq(x_range[1], x_range[2], length.out = 100)
      robust_line_data <- data.frame(logFC_A = x_seq, logFC_B_predicted = predict(rlm_fit, newdata = data.frame(logFC_A = x_seq)))
    }, error = function(e) {
      warning(paste("RLM fitting failed for", gene_highlight, labA, "vs", labB, ":", e$message))
      robust_line_data <<- NULL
    })
  }
  
  pal <- setNames(c(colBoth, colA, colB, colNS), c("Both", labA, labB, "Not significant"))
  
  p <- ggplot(dge_tab, aes(x = logFC_A, y = logFC_B)) +
    geom_point(data = . %>% filter(significant_cat == "Not significant"), aes(colour = significant_cat), alpha = .5, size = 1.5) +
    geom_point(data = . %>% filter(significant_cat != "Not significant"), aes(colour = significant_cat), alpha = .7, size = 2) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = .3, color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = .3, color = "grey50") +
    ggpubr::stat_cor(aes(group = 1), method = "pearson", label.x.npc = "right", label.y.npc = "top", size = 3, color = "black") +
    ggpubr::stat_cor(aes(group = 1), method = "spearman", label.x.npc = "right", label.y.npc = 0.93, size = 3, color = "darkblue") +
    geom_point(data = dge_tab %>% filter(gene_symbol == gene_highlight), inherit.aes = FALSE, aes(x = logFC_A, y = logFC_B), colour = "black", fill = "gold", shape = 21, size = 3.5, stroke = 1) +
    ggrepel::geom_text_repel(data  = dge_tab %>% filter(gene_symbol == gene_highlight), inherit.aes = FALSE, aes(label = gene_symbol, x = logFC_A, y = logFC_B), colour = "black", fontface = "bold", size = 3, arrow = grid::arrow(length = unit(0.02, "npc"), type = "closed"), max.overlaps = Inf, box.padding = 0.4, point.padding = 0.15) +
    labs(title = paste0(gene_highlight, ": ", labA, " vs ", labB), x = paste0("logFC (", labA, ")"), y = paste0("logFC (", labB, ")")) +
    scale_colour_manual(values = pal, breaks = names(pal), name = "Significance", drop = FALSE) +
    theme_bw(base_size=10) +
    theme(plot.title = element_text(hjust = .5, face = "bold"), legend.position = "bottom", legend.title = element_text(face="bold"))
  
  if (!is.null(robust_line_data)) {
    p <- p + geom_line(data = robust_line_data, aes(x = logFC_A, y = logFC_B_predicted), colour = "firebrick", linewidth = 0.7, linetype = "solid")
  } else {
    p <- p + geom_smooth(aes(group=1), method="lm", se=FALSE, colour="black", linewidth=0.4, linetype="dotted", formula=y~x)
  }
  return(p)
}

corr_plot_combos <- list(
  list(df1_key = "Fas_KO", df2_key = "Fas_KD", highlight_gene = "FAS", labA = "KO", labB = "KD"),
  list(df1_key = "RC3H1_KO", df2_key = "RC3H1_KD", highlight_gene = "RC3H1", labA = "KO", labB = "KD"),
  list(df1_key = "PTPN2_KO", df2_key = "PTPN2_KD", highlight_gene = "PTPN2", labA = "KO", labB = "KD"),
  list(df1_key = "SUV39H1_KO", df2_key = "SUV39H1_KD", highlight_gene = "SUV39H1", labA = "KO", labB = "KD"),
  list(df1_key = "MED12_KO", df2_key = "MED12_KD", highlight_gene = "MED12", labA = "KO", labB = "KD"),
  list(df1_key = "RASA2_KO", df2_key = "RASA2_KD", highlight_gene = "RASA2", labA = "KO", labB = "KD")
)

correlation_plots_list <- lapply(corr_plot_combos, function(p) {
  if (p$df1_key %in% names(all_results_topTable) && p$df2_key %in% names(all_results_topTable)) {
    create_corr_plot_robust_rlm(
      df1 = all_results_topTable[[p$df1_key]],
      df2 = all_results_topTable[[p$df2_key]],
      gene_highlight = p$highlight_gene,
      labA = p$labA,
      labB = p$labB
    )
  } else {
    warning(paste("Missing results for correlation plot:", p$df1_key, "or", p$df2_key))
    NULL
  }
})
correlation_plots_list <- Filter(Negate(is.null), correlation_plots_list)

if(length(correlation_plots_list) > 0) {
  combined_corr_plot <- wrap_plots(correlation_plots_list, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path("plots_dge", paste0(EXPT_NAME, "_combined_correlation_plots_KOvKD.png")), combined_corr_plot,
         width = 10, height = 4 + 5 * ceiling(length(correlation_plots_list)/2))
}

# Heatmap for Double Differential Results
create_double_diff_heatmap <- function(all_results, targets) {
  sig_genes_union <- c()
  
  for (target in targets) {
    dd_contrast <- paste0(target, "_KDvsKO")
    if (dd_contrast %in% names(all_results)) {
      sig_genes <- all_results[[dd_contrast]] %>%
        filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) %>%
        pull(gene_symbol) %>%
        na.omit()
      sig_genes_union <- union(sig_genes_union, sig_genes)
    }
  }
  
  if (length(sig_genes_union) > 0) {
    mat <- matrix(NA, nrow = length(sig_genes_union), ncol = length(targets))
    rownames(mat) <- sig_genes_union
    colnames(mat) <- targets
    
    for (i in seq_along(targets)) {
      dd_contrast <- paste0(targets[i], "_KDvsKO")
      if (dd_contrast %in% names(all_results)) {
        result_df <- all_results[[dd_contrast]]
        matched_idx <- match(sig_genes_union, result_df$gene_symbol)
        mat[, i] <- result_df$logFC[matched_idx]
      }
    }
    
    mat[is.na(mat)] <- 0
    
    if (nrow(mat) > 50) {
      mat <- mat[order(rowSums(abs(mat)), decreasing = TRUE)[1:50], ]
    }
    
    hm <- Heatmap(mat,
                  name = "logFC\n(KD-KO)",
                  cluster_rows = TRUE,
                  cluster_columns = FALSE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 10),
                  column_title = "Double Differential: KD vs KO",
                  column_title_gp = gpar(fontsize = 12, fontface = "bold"))
    
    png(file.path("plots_dge", paste0(EXPT_NAME, "_double_differential_heatmap.png")), 
        width = 8, height = 10, units = "in", res = 300)
    draw(hm)
    dev.off()
  }
}

create_double_diff_heatmap(all_results_topTable, targets)

# Summary Statistics
summary_stats <- data.frame(
  Contrast = character(),
  Total_Tested = integer(),
  Sig_Up = integer(),
  Sig_Down = integer(),
  Total_Sig = integer(),
  stringsAsFactors = FALSE
)

for (contrast_name in names(all_results_topTable)) {
  result <- all_results_topTable[[contrast_name]]
  sig_up <- sum(result$adj.P.Val < P_VALUE_CUTOFF & result$logFC > LOGFC_CUTOFF, na.rm = TRUE)
  sig_down <- sum(result$adj.P.Val < P_VALUE_CUTOFF & result$logFC < -LOGFC_CUTOFF, na.rm = TRUE)
  
  summary_stats <- rbind(summary_stats, data.frame(
    Contrast = contrast_name,
    Total_Tested = nrow(result),
    Sig_Up = sig_up,
    Sig_Down = sig_down,
    Total_Sig = sig_up + sig_down
  ))
}

if (!dir.exists("results_dge")) dir.create("results_dge", recursive = TRUE)
write_csv(summary_stats, file.path("results_dge", paste0(EXPT_NAME, "_summary_statistics.csv")))

# Venn Diagrams
get_sig <- function(tbl) {
  if (is.null(tbl)) return(character(0))
  tbl %>%
    filter(adj.P.Val < P_VALUE_CUTOFF, abs(logFC) > LOGFC_CUTOFF) %>%
    pull(gene_symbol) %>%
    na.omit() %>%
    unique()
}

for (tg in targets) {
  ko_name <- paste0(tg, "_KO")
  kd_name <- paste0(tg, "_KD")
  
  ko_sig <- get_sig(all_results_topTable[[ko_name]])
  kd_sig <- get_sig(all_results_topTable[[kd_name]])
  
  ko_sig <- unique(c(ko_sig, tg))
  
  venn_sets <- list(KO = ko_sig, KD = kd_sig)
  
  venn_plot <- ggVennDiagram(venn_sets, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste("KO vs KD overlap –", tg)) +
    theme(legend.position = "none")
  
  ggsave(file.path("plots_dge", paste0(EXPT_NAME, "_", tg, "_KOvsKD_venn.png")),
         venn_plot, width = 5, height = 4)
}

# Final Summary
cat("\n========================================\n")
cat("Analysis for experiment", EXPT_NAME, "complete!\n")
cat("========================================\n\n")
cat("Output directories:\n")
cat("- Results CSVs: 'results_dge/'\n")
cat("- GO/KEGG/MSigDB results: 'results_go/'\n")
cat("- Enrichment comparisons: 'results_enrichment_comparison/'\n")
cat("- All plots: 'plots_dge/'\n")
cat("\nKey analyses performed:\n")
cat("- Standard differential expression (KO vs NTC, KD vs NTC)\n")
cat("- Double differential analysis (KD vs KO)\n")
cat("- GO, KEGG, and MSigDB Hallmark enrichment analysis\n")
cat("- Cross-contrast enrichment comparisons\n")
cat("- Identification of commonly enriched pathways\n")
cat("- Comprehensive visualizations (volcano, correlation, heatmaps, venn)\n")
cat("\nConfiguration used:\n")
cat("- P-value cutoff:", P_VALUE_CUTOFF, "\n")
cat("- LogFC cutoff:", LOGFC_CUTOFF, "\n")
cat("- Min genes for enrichment:", MIN_GENES_FOR_ENRICHMENT, "\n")
cat("\n========================================\n")