library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Configuration
EXPERIMENT_NAME <- "Co065"
proximal_distance_kb <- 100

# Create output directories
if (!dir.exists("Co065_offtarget_results")) dir.create("Co065_offtarget_results", recursive = TRUE)
if (!dir.exists("Co065_heatmap_results")) dir.create("Co065_heatmap_results", recursive = TRUE)

cat("CO065 UNIFIED OFF-TARGET ANALYSIS\n\n")

# MAPPING CO065 DGE FILES TO GUIDES
guide_to_dge_mapping <- list(
  # FAS files  
  "FAS_KD" = "DE_topTable_Co065_Fas_KD.csv",
  "FAS_KO" = "DE_topTable_Co065_Fas_KO.csv",
  
  # MED12 files
  "MED12_KD" = "DE_topTable_Co065_MED12_KD.csv", 
  "MED12_KO" = "DE_topTable_Co065_MED12_KO.csv",
  
  # PTPN2 files
  "PTPN2_KD" = "DE_topTable_Co065_PTPN2_KD.csv",
  "PTPN2_KO" = "DE_topTable_Co065_PTPN2_KO.csv",
  
  # RASA2 files
  "RASA2_KD" = "DE_topTable_Co065_RASA2_KD.csv",
  "RASA2_KO" = "DE_topTable_Co065_RASA2_KO.csv",
  
  # RC3H1 files
  "RC3H1_KD" = "DE_topTable_Co065_RC3H1_KD.csv",
  "RC3H1_KO" = "DE_topTable_Co065_RC3H1_KO.csv",
  
  # SUV39H1 files
  "SUV39H1_KD" = "DE_topTable_Co065_SUV39H1_KD.csv",
  "SUV39H1_KO" = "DE_topTable_Co065_SUV39H1_KO.csv"
)

# PROMOTER ANNOTATION SETUP
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
proms <- promoters(genes, upstream = 1000, downstream = 1000)

if (length(proms) > 0) {
  entrez_ids <- names(proms)
  mcols(proms)$entrez_id <- entrez_ids
  mcols(proms)$symbol <- mapIds(
    org.Hs.eg.db, keys = entrez_ids,
    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
  )
}

# HELPER FUNCTIONS
parse_locus <- function(locus_vec) {
  tibble(
    chr    = str_extract(locus_vec, "^chr[^:]+"),
    strand = case_when(str_detect(locus_vec, ":[+]") ~ "+",
                       str_detect(locus_vec, ":-")  ~ "-",
                       TRUE                         ~ "*"),
    pos    = as.integer(str_remove(locus_vec, "^chr[^:]+:[+-]"))
  )
}

get_ontarget_location <- function(df) {
  if (nrow(df) > 0) {
    return(df$Locus[1])
  } else {
    return(NA)
  }
}

find_proximal_genes <- function(ontarget_locus, distance_kb = 100) {
  if (is.na(ontarget_locus)) {
    return(tibble(gene_symbol = character(), distance = numeric(), gene_type = character()))
  }
  
  parsed_locus <- parse_locus(ontarget_locus)
  if (nrow(parsed_locus) == 0 || is.na(parsed_locus$pos[1])) {
    return(tibble(gene_symbol = character(), distance = numeric(), gene_type = character()))
  }
  
  # Create GRanges object for on-target location
  ontarget_gr <- GRanges(
    seqnames = parsed_locus$chr[1],
    ranges = IRanges(start = parsed_locus$pos[1], width = 1),
    strand = parsed_locus$strand[1]
  )
  
  # Find genes within distance
  distance_bp <- distance_kb * 1000
  expanded_region <- resize(ontarget_gr, width = distance_bp * 2, fix = "center")
  
  overlaps <- findOverlaps(expanded_region, genes)
  
  if (length(overlaps) > 0) {
    overlapping_genes <- genes[subjectHits(overlaps)]
    
    # Calculate distances from on-target site to gene start positions
    gene_distances <- distance(ontarget_gr, overlapping_genes)
    
    # Get gene symbols
    gene_symbols <- mapIds(
      org.Hs.eg.db, 
      keys = names(overlapping_genes),
      column = "SYMBOL", 
      keytype = "ENTREZID", 
      multiVals = "first"
    )
    
    result <- tibble(
      gene_symbol = as.character(gene_symbols),
      distance = as.numeric(gene_distances),
      gene_type = "proximal_gene"
    ) %>%
      filter(!is.na(gene_symbol)) %>%
      arrange(distance)
    
    return(result)
  } else {
    return(tibble(gene_symbol = character(), distance = numeric(), gene_type = character()))
  }
}

# PROCESS OFF-TARGET DATA
cat("Processing off-target data files...\n")

offtarget_files <- list.files("offtarget_data", pattern = "\\.csv$", full.names = TRUE)

process_offtarget_file <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  target_name <- str_remove(basename(file_path), "\\.csv$")
  
  if (nrow(df) == 0) {
    cat("  ", target_name, ": No data\n")
    return(NULL)
  }
  
  # Get on-target location for proximal analysis
  ontarget_locus <- get_ontarget_location(df)
  
  # Process off-target sites
  locus_parsed <- parse_locus(df$Locus)
  df_with_coords <- bind_cols(df, locus_parsed)
  
  offtarget_gr <- GRanges(
    seqnames = df_with_coords$chr,
    ranges = IRanges(start = df_with_coords$pos, width = 1),
    strand = df_with_coords$strand
  )
  
  # Determine if this is KD (CRISPRoff) or KO (Cas9) experiment
  is_ko_experiment <- str_detect(target_name, regex("KO|maya|cas9", ignore_case = TRUE))
  
  if (is_ko_experiment) {
    # For KO experiments: use Gene column from IDT off-target file
    # Distance from promoter doesn't matter - only if edit is on an exon (IDT handles this)
    offtarget_with_genes <- df %>%
      filter(!is.na(Gene) & Gene != "" & Gene != " ")
    
    if (nrow(offtarget_with_genes) > 0) {
      offtarget_results <- tibble(
        target_file = target_name,
        locus = offtarget_with_genes$Locus,
        gene_symbol = offtarget_with_genes$Gene,
        entrez_id = NA_character_,
        mismatch = as.character(offtarget_with_genes$Score),
        gene_type = "Gene",
        analysis_type = "offtarget",
        distance = NA_real_
      ) %>%
        filter(!is.na(gene_symbol) & gene_symbol != "")
    } else {
      offtarget_results <- tibble(
        target_file = character(),
        locus = character(),
        gene_symbol = character(),
        entrez_id = character(),
        mismatch = character(),
        gene_type = character(),
        analysis_type = character(),
        distance = numeric()
      )
    }
  } else {
    # For KD experiments: use promoter overlap logic
    # Distance from promoter matters for CRISPRoff knockdown
    overlaps <- findOverlaps(offtarget_gr, proms)
    
    if (length(overlaps) > 0) {
      df_hits <- df_with_coords[queryHits(overlaps), ]
      prom_hits <- proms[subjectHits(overlaps)]
      
      offtarget_results <- tibble(
        target_file = target_name,
        locus = df_hits$Locus,
        gene_symbol = as.character(mcols(prom_hits)$symbol),
        entrez_id = as.character(mcols(prom_hits)$entrez_id),
        mismatch = as.character(df_hits$Score),
        gene_type = "Promoter",
        analysis_type = "offtarget",
        distance = NA_real_
      ) %>%
        filter(!is.na(gene_symbol))
    } else {
      offtarget_results <- tibble(
        target_file = character(),
        locus = character(),
        gene_symbol = character(),
        entrez_id = character(),
        mismatch = character(),
        gene_type = character(),
        analysis_type = character(),
        distance = numeric()
      )
    }
  }
  
  # Find proximal genes
  proximal_results <- find_proximal_genes(ontarget_locus, proximal_distance_kb)
  
  if (nrow(proximal_results) > 0) {
    proximal_results <- proximal_results %>%
      mutate(
        target_file = target_name,
        locus = ontarget_locus,
        entrez_id = NA_character_,
        mismatch = as.character(distance),  # Store distance as mismatch for sorting
        analysis_type = "proximal"
      ) %>%
      dplyr::select(target_file, locus, gene_symbol, entrez_id, mismatch, gene_type, analysis_type, distance)
  }
  
  # Combine results
  combined_results <- bind_rows(offtarget_results, proximal_results)
  
  if (nrow(combined_results) > 0) {
    # Save individual file results
    output_file <- file.path("Co065_offtarget_results", paste0(target_name, "_annotated.csv"))
    write_csv(combined_results, output_file)
    
    cat("  ", target_name, ": ", nrow(offtarget_results), " off-target genes, ", 
        nrow(proximal_results), " proximal genes\n")
    
    return(combined_results)
  } else {
    cat("  ", target_name, ": No genes found\n")
    return(NULL)
  }
}

# Process all files
all_offtarget_data <- map_dfr(offtarget_files, process_offtarget_file)

cat("Off-target annotation complete!\n\n")

# READ DGE RESULTS
cat("Reading DGE results...\n")

read_dge_results <- function() {
  dge_files <- list.files("results_Co065", pattern = "^DE_topTable_Co065_.*\\.(csv|CSV)$", full.names = TRUE)
  
  dge_data <- map_dfr(dge_files, function(file) {
    df <- read_csv(file, show_col_types = FALSE)
    df$dge_file <- basename(file)
    return(df)
  })
  
  return(dge_data)
}

dge_data <- read_dge_results()

# CREATE GUIDE-SPECIFIC DATA MATRICES
cat("Creating guide-specific data matrices...\n")

# Create individual datasets for each guide
create_guide_specific_data <- function(guide_mapping, offtarget_data, dge_data) {
  guide_matrices <- list()
  
  # Create mapping between guide names and offtarget files
  # For pooled experiments: combine offtarget data from all available guides for that target
  # For individual guides: use specific guide data when available
  
  guide_to_offtarget_mapping <- list(
    # FAS: KD uses CRISPRoff guides, KO uses separate KO file
    "FAS_KD" = c("FAS_CRISPRoff_guide1", "FAS_CRISPRoff_guide2", "FAS_CRISPRoff_guide3"),
    "FAS_KO" = c("FAS_KO_g1"),
    
    # PTPN2: KD uses CRISPRoff guides, KO uses separate KO file
    "PTPN2_KD" = c("PTPN2_CRISPRoff_guide3", "PTPN2_CRISPRoff_guide5", "PTPN2_CRISPRoff_guide6"),
    "PTPN2_KO" = c("PTPN2_KO_g1"),
    
    # RASA2: KD uses CRISPRoff guides, KO uses separate KO file
    "RASA2_KD" = c("RASA2_CRISPRoff_guide1", "RASA2_CRISPRoff_guide2", "RASA2_CRISPRoff_guide4"),
    "RASA2_KO" = c("RASA2_KO_g1"),
    
    # RC3H1: KD uses CRISPRoff guides, KO uses separate KO file
    "RC3H1_KD" = c("RC3H1_CRISPRoff_guide1", "RC3H1_CRISPRoff_guide2", "RC3H1_CRISPRoff_guide3"),
    "RC3H1_KO" = c("RC3H1_KO_g1"),
    
    # SUV39H1: KD uses CRISPRoff guides, KO uses separate KO file
    "SUV39H1_KD" = c("SUV39H1_CRISPRoff_guide2", "SUV39H1_CRISPRoff_guide3"),
    "SUV39H1_KO" = c("SUV39H1_KO_g2"),
    
    # MED12: KD uses CRISPRoff guides, KO uses separate KO file
    "MED12_KD" = c("MED12_CRISPRoff_guide1", "MED12_CRISPRoff_guide2", "MED12_CRISPRofff_guide3"),
    "MED12_KO" = c("MED12_KO_maya")
  )
  
  for (guide_name in names(guide_mapping)) {
    dge_file <- guide_mapping[[guide_name]]
    
    # Extract target gene from guide name
    target_gene <- str_extract(guide_name, "^[^_]+")
    target_gene <- str_to_upper(target_gene)
    
    # Find corresponding offtarget files using mapping
    offtarget_files <- guide_to_offtarget_mapping[[guide_name]]
    
    if (is.null(offtarget_files)) {
      cat("No offtarget mapping found for guide:", guide_name, "\n")
      next
    }
    
    # Get genes for this guide from offtarget data (combining multiple files for pooled experiments)
    guide_genes_raw <- offtarget_data %>%
      filter(target_file %in% offtarget_files) %>%
      dplyr::select(gene_symbol, analysis_type, gene_type, mismatch, distance) %>%
      distinct()
    
    # Add the target gene to the analysis
    target_gene_entry <- tibble(
      gene_symbol = target_gene,
      analysis_type = "target",
      gene_type = "Target",
      mismatch = "0",
      distance = 0
    )
    
    guide_genes_raw <- bind_rows(target_gene_entry, guide_genes_raw)
    
    # For pooled experiments with multiple offtarget files, we need to handle duplicates
    # Keep the best (lowest) mismatch score for offtarget genes
    if (length(offtarget_files) > 1) {
      guide_genes_raw <- guide_genes_raw %>%
        group_by(gene_symbol, analysis_type) %>%
        arrange(as.numeric(ifelse(mismatch == "N/A" | is.na(mismatch), 999, mismatch))) %>%
        slice_head(n = 1) %>%
        ungroup()
    }
    
    # Identify genes that appear in both categories
    both_genes <- guide_genes_raw %>%
      group_by(gene_symbol) %>%
      summarise(is_both = n_distinct(analysis_type) > 1, .groups = "drop") %>%
      filter(is_both) %>%
      pull(gene_symbol)
    
    # Process genes: keep all genes but mark ones that are "both"
    guide_genes <- guide_genes_raw %>%
      mutate(
        # Mark genes that appear in both categories
        analysis_type = ifelse(gene_symbol %in% both_genes, "both", analysis_type),
        mismatch_numeric = ifelse(mismatch == "N/A" | is.na(mismatch), 0, 
                                 ifelse(is.na(as.numeric(mismatch)), 999, as.numeric(mismatch))),
        distance_numeric = ifelse(is.na(distance), 999, distance),
        category_order = case_when(
          analysis_type == "target" ~ 0,
          analysis_type == "offtarget" ~ 1,
          analysis_type == "both" ~ 2,
          analysis_type == "proximal" ~ 3,
          TRUE ~ 4
        ),
        sort_value = case_when(
          analysis_type == "target" ~ 0,
          analysis_type == "offtarget" | analysis_type == "both" ~ mismatch_numeric,
          analysis_type == "proximal" ~ distance_numeric,
          TRUE ~ 999
        )
      ) %>%
      # For genes in both categories, keep the offtarget version for metadata
      group_by(gene_symbol) %>%
      arrange(analysis_type == "proximal") %>%  # offtarget/both comes first
      slice_head(n = 1) %>%
      ungroup()
    
    if (nrow(guide_genes) == 0) {
      cat("No genes found for guide:", guide_name, "(offtarget files:", paste(offtarget_files, collapse = ", "), ")\n")
      next
    }
    
    # Check if we have genes
    if (nrow(guide_genes) == 0) {
      cat("No genes found for guide:", guide_name, "\n")
      next
    }
    
    # Get DGE data for this guide
    guide_dge <- dge_data %>%
      filter(dge_file == !!dge_file) %>%
      filter(!is.na(gene_symbol)) %>%
      filter(gene_symbol %in% guide_genes$gene_symbol)
    
    if (nrow(guide_dge) == 0) {
      cat("No DGE data found for guide:", guide_name, "- will include all genes with 0 logFC\n")
      # Create empty DGE data frame but continue with analysis
      guide_dge <- tibble(gene_symbol = character(), logFC = numeric(), adj.P.Val = numeric())
    }
    
    # Determine guide type
    if (str_detect(guide_name, "pool|single")) {
      guide_type <- "KD"
    } else if (str_detect(guide_name, "KO|maya")) {
      guide_type <- "KO"
    } else {
      guide_type <- "KD"  # Default to KD for other guides
    }
    
    # Create single-row matrix (guide x genes) - handle duplicates in DGE data too
    guide_dge_unique <- guide_dge %>%
      group_by(gene_symbol) %>%
      summarise(logFC = dplyr::first(logFC), adj.P.Val = dplyr::first(adj.P.Val), .groups = "drop")
    
    lfc_vector <- guide_dge_unique$logFC
    names(lfc_vector) <- guide_dge_unique$gene_symbol
    
    # Create custom gene ordering: offtarget (by score) -> proximal (by distance)
    gene_order <- guide_genes %>%
      arrange(category_order, sort_value) %>%
      pull(gene_symbol)
    
    # Ensure all genes from offtarget analysis are included
    all_genes <- gene_order
    missing_genes <- setdiff(all_genes, names(lfc_vector))
    
    # Add missing genes with 0 logFC
    if (length(missing_genes) > 0) {
      missing_values <- rep(0, length(missing_genes))
      names(missing_values) <- missing_genes
      lfc_vector <- c(lfc_vector, missing_values)
    }
    
    # Create matrix with single row using custom ordering
    lfc_matrix <- matrix(lfc_vector[all_genes], nrow = 1, dimnames = list(guide_name, all_genes))
    
    # Create significance vector
    sig_vector <- guide_dge_unique$adj.P.Val < 0.05
    names(sig_vector) <- guide_dge_unique$gene_symbol
    
    # Add missing genes with FALSE significance
    if (length(missing_genes) > 0) {
      missing_sig <- rep(FALSE, length(missing_genes))
      names(missing_sig) <- missing_genes
      sig_vector <- c(sig_vector, missing_sig)
    }
    
    # Create significance matrix
    sig_matrix <- matrix(sig_vector[all_genes], nrow = 1, dimnames = list(guide_name, all_genes))
    
    # Create gene annotation with ordered genes
    gene_annotation <- guide_genes %>%
      arrange(match(gene_symbol, all_genes)) %>%
      dplyr::select(gene_symbol, analysis_type, gene_type, mismatch_numeric, distance_numeric)
    
    guide_matrices[[guide_name]] <- list(
      guide_name = guide_name,
      target_gene = target_gene,
      guide_type = guide_type,
      lfc_matrix = lfc_matrix,
      sig_matrix = sig_matrix,
      gene_annotation = gene_annotation,
      offtarget_files = paste(offtarget_files, collapse = ", ")
    )
    
    cat("Processed guide:", guide_name, "- Genes:", ncol(lfc_matrix), "\n")
  }
  
  return(guide_matrices)
}

guide_matrices <- create_guide_specific_data(guide_to_dge_mapping, all_offtarget_data, dge_data)

# CREATE HEATMAPS
cat("Generating individual guide heatmaps...\n")

create_guide_heatmaps <- function(guide_matrices, output_dir = "Co065_heatmap_results") {
  
  heatmap_list <- list()
  
  for (guide_name in names(guide_matrices)) {
    guide_data <- guide_matrices[[guide_name]]
    
    lfc_matrix <- guide_data$lfc_matrix
    sig_matrix <- guide_data$sig_matrix
    gene_annotation <- guide_data$gene_annotation
    target_gene <- guide_data$target_gene
    guide_type <- guide_data$guide_type
    
    if (ncol(lfc_matrix) > 0 && nrow(lfc_matrix) > 0) {
      
      # Create significance annotation
      sig_text <- sig_matrix
      sig_text[sig_matrix] <- "*"
      sig_text[!sig_matrix] <- ""
      
      # Transpose matrix so genes are on x-axis
      lfc_matrix_t <- t(lfc_matrix)
      sig_text_t <- t(sig_text)
      
      # Create heatmap - genes on x-axis
      ht <- Heatmap(
        lfc_matrix_t,
        name = "logFC",
        col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
        
        # Row (gene) settings
        row_title = "Genes",
        row_title_side = "left",
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        show_row_names = TRUE,
        
        # Column (guide) settings - single column
        column_title = paste0(guide_name, " (", guide_type, ")"),
        column_title_side = "top",
        column_names_rot = 0,
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        show_column_names = TRUE,
        
        # Add row annotation for gene types
        left_annotation = rowAnnotation(
          Analysis = gene_annotation$analysis_type,
          col = list(Analysis = c("target" = "red", "offtarget" = "orange", "proximal" = "green", "both" = "purple")),
          annotation_name_side = "bottom",
          width = unit(0.5, "cm")
        ),
        
        # Add significance stars
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (sig_text_t[i, j] == "*") {
            grid.text("*", x, y, gp = gpar(fontsize = 8, col = "black"))
          }
        },
        
        # No clustering - use custom ordering
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        # Size calculations
        width = unit(4, "cm"),  # Fixed width for single column
        height = unit(max(10, nrow(lfc_matrix_t) * 0.3), "cm"),
        
        # Add cell borders
        rect_gp = gpar(col = "white", lwd = 0.5)
      )
      
      # Save heatmap
      safe_guide_name <- str_replace_all(guide_name, "[^A-Za-z0-9_-]", "_")
      output_file <- file.path(output_dir, paste0("Co065_", safe_guide_name, "_heatmap.pdf"))
      output_file_png <- file.path(output_dir, paste0("Co065_", safe_guide_name, "_heatmap.png"))
      
      # Calculate dynamic figure dimensions
      fig_width <- 8  # Fixed width for single column
      fig_height <- max(12, nrow(lfc_matrix_t) * 0.25 + 4)
      
      pdf(output_file, width = fig_width, height = fig_height)
      draw(ht)
      dev.off()
      
      png(output_file_png, width = fig_width * 100, height = fig_height * 100, res = 100)
      draw(ht)
      dev.off()
      
      heatmap_list[[guide_name]] <- ht
      
      cat("Created heatmap for", guide_name, "- saved to", output_file, "and", output_file_png, "\n")
      cat("  Target:", target_gene, "| Type:", guide_type, "| Genes:", ncol(lfc_matrix), "\n")
    }
  }
  
  return(heatmap_list)
}

# Export heatmap data
export_guide_data <- function(guide_matrices, output_dir = "Co065_heatmap_results") {
  
  for (guide_name in names(guide_matrices)) {
    guide_data <- guide_matrices[[guide_name]]
    
    lfc_matrix <- guide_data$lfc_matrix
    gene_annotation <- guide_data$gene_annotation
    target_gene <- guide_data$target_gene
    guide_type <- guide_data$guide_type
    
    # Export logFC data for this guide
    lfc_data <- data.frame(
      gene_symbol = colnames(lfc_matrix),
      logFC = as.vector(lfc_matrix[1, ]),
      significance = as.vector(guide_data$sig_matrix[1, ]),
      stringsAsFactors = FALSE
    ) %>%
      left_join(gene_annotation %>% dplyr::select(gene_symbol, analysis_type, gene_type), by = "gene_symbol")
    
    safe_guide_name <- str_replace_all(guide_name, "[^A-Za-z0-9_-]", "_")
    lfc_file <- file.path(output_dir, paste0("Co065_", safe_guide_name, "_logFC_data.csv"))
    write_csv(lfc_data, lfc_file)
    
    # Export guide info
    guide_info <- tibble(
      guide_name = guide_name,
      target_gene = target_gene,
      guide_type = guide_type,
      offtarget_files = guide_data$offtarget_files,
      total_genes = ncol(lfc_matrix),
      target_genes = sum(gene_annotation$analysis_type == "target"),
      offtarget_genes = sum(gene_annotation$analysis_type == "offtarget"),
      proximal_genes = sum(gene_annotation$analysis_type == "proximal"),
      both_genes = sum(gene_annotation$analysis_type == "both")
    )
    
    guide_file <- file.path(output_dir, paste0("Co065_", safe_guide_name, "_guide_info.csv"))
    write_csv(guide_info, guide_file)
    
    cat("Data for", guide_name, "saved to", lfc_file, "and", guide_file, "\n")
  }
}

# Create summary tables
create_guide_summaries <- function(guide_matrices, output_dir = "Co065_heatmap_results") {
  
  summary_list <- list()
  all_summaries <- list()
  
  for (guide_name in names(guide_matrices)) {
    guide_data <- guide_matrices[[guide_name]]
    gene_annotation <- guide_data$gene_annotation
    
    # Create summary of genes by analysis type for this guide
    gene_summary <- gene_annotation %>%
      count(analysis_type, gene_type, sort = TRUE) %>%
      mutate(guide_name = guide_name, target_gene = guide_data$target_gene)
    
    safe_guide_name <- str_replace_all(guide_name, "[^A-Za-z0-9_-]", "_")
    summary_file <- file.path(output_dir, paste0("Co065_", safe_guide_name, "_gene_summary.csv"))
    write_csv(gene_summary, summary_file)
    
    summary_list[[guide_name]] <- gene_summary
    all_summaries[[guide_name]] <- gene_summary
    
    cat("Gene summary for", guide_name, "saved to", summary_file, "\n")
  }
  
  # Create overall summary across all guides
  overall_summary <- bind_rows(all_summaries) %>%
    group_by(target_gene, analysis_type, gene_type) %>%
    summarise(total_guides = n_distinct(guide_name), 
              total_genes = sum(n), 
              .groups = "drop") %>%
    arrange(target_gene, analysis_type)
  
  overall_file <- file.path(output_dir, "Co065_all_guides_summary.csv")
  write_csv(overall_summary, overall_file)
  cat("Overall summary saved to", overall_file, "\n")
  
  return(list(individual = summary_list, overall = overall_summary))
}

# RUN ANALYSIS
heatmaps <- create_guide_heatmaps(guide_matrices)

cat("\nExporting guide data...\n")
export_guide_data(guide_matrices)

cat("Creating summary tables...\n")
summaries <- create_guide_summaries(guide_matrices)

cat("\nCO065 ANALYSIS COMPLETE\n")
cat("Generated", length(heatmaps), "individual guide heatmaps:\n")
for (guide_name in names(heatmaps)) {
  guide_data <- guide_matrices[[guide_name]]
  cat(sprintf("  %-25s: %s (%s) - %d genes\n", 
              guide_name, 
              guide_data$target_gene,
              guide_data$guide_type,
              ncol(guide_data$lfc_matrix)))
}
cat("All results saved to Co065_offtarget_results/ and Co065_heatmap_results/ directories\n")

# Print summary by target gene
cat("\nSUMMARY BY TARGET GENE\n")
target_summary <- map_dfr(guide_matrices, function(x) {
  tibble(
    target_gene = x$target_gene,
    guide_name = x$guide_name,
    guide_type = x$guide_type,
    total_genes = ncol(x$lfc_matrix),
    offtarget_genes = sum(x$gene_annotation$analysis_type == "offtarget"),
    proximal_genes = sum(x$gene_annotation$analysis_type == "proximal")
  )
}) %>%
  arrange(target_gene, guide_name)

for (target in unique(target_summary$target_gene)) {
  target_guides <- target_summary %>% filter(target_gene == target)
  cat(sprintf("\n%s:\n", target))
  for (i in 1:nrow(target_guides)) {
    cat(sprintf("  %-20s (%s): %d total genes (%d offtarget, %d proximal)\n",
                target_guides$guide_name[i],
                target_guides$guide_type[i], 
                target_guides$total_genes[i],
                target_guides$offtarget_genes[i],
                target_guides$proximal_genes[i]))
  }
}