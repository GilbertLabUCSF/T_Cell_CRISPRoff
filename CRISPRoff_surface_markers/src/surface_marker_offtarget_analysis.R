# Surface Marker Off-target Analysis with Cross-reference and Heatmaps
library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Configuration
EXPERIMENT_NAME <- "Surface_Markers"
proximal_distance_kb <- 100

# Create output directories
if (!dir.exists("offtarget_results")) dir.create("offtarget_results", recursive = TRUE)
if (!dir.exists("heatmap_results")) dir.create("heatmap_results", recursive = TRUE)

cat("SURFACE MARKER OFF-TARGET ANALYSIS\n\n")

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

# HELPER FUNCTIONS (same as Co062)
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
  
  # Define search window
  distance_bp <- distance_kb * 1000
  
  # Find overlapping promoters within the window
  window_start <- start(ontarget_gr) - distance_bp
  window_end <- start(ontarget_gr) + distance_bp
  
  # Create window GRanges
  window_gr <- GRanges(
    seqnames = seqnames(ontarget_gr),
    ranges = IRanges(start = window_start, end = window_end),
    strand = "*"
  )
  
  # Find overlapping promoters
  overlaps <- findOverlaps(window_gr, proms)
  
  if (length(overlaps) > 0) {
    promoter_hits <- proms[subjectHits(overlaps)]
    
    # Calculate distances
    distances <- abs(start(promoter_hits) - start(ontarget_gr))
    
    proximal_genes <- tibble(
      gene_symbol = as.character(mcols(promoter_hits)$symbol),
      distance = distances,
      gene_type = "Promoter"
    ) %>%
      filter(!is.na(gene_symbol)) %>%
      arrange(distance)
  } else {
    proximal_genes <- tibble(gene_symbol = character(), distance = numeric(), gene_type = character())
  }
  
  return(proximal_genes)
}

# PROCESS OFF-TARGET FILES
process_offtarget_file <- function(file_path, target_name, target_gene_symbol) {
  cat("\nProcessing:", target_name, "\n")
  
  # Read off-target file
  df <- read_csv(file_path, show_col_types = FALSE)
  
  if (nrow(df) == 0) {
    cat("No off-targets found in file.\n")
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
  
  # For CRISPRoff: use promoter overlap logic
  overlaps <- findOverlaps(offtarget_gr, proms)
  
  if (length(overlaps) > 0) {
    df_hits <- df_with_coords[queryHits(overlaps), ]
    promoter_hits <- proms[subjectHits(overlaps)]
    
    offtarget_results <- tibble(
      target_file = target_name,
      locus = df_hits$Locus,
      gene_symbol = as.character(mcols(promoter_hits)$symbol),
      entrez_id = as.character(mcols(promoter_hits)$entrez_id),
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
  
  # Find proximal genes
  proximal_results <- find_proximal_genes(ontarget_locus, proximal_distance_kb)
  
  if (nrow(proximal_results) > 0) {
    proximal_results <- proximal_results %>%
      mutate(
        target_file = target_name,
        locus = ontarget_locus,
        entrez_id = NA_character_,
        mismatch = as.character(distance),
        analysis_type = "proximal"
      ) %>%
      dplyr::select(target_file, locus, gene_symbol, entrez_id, mismatch, gene_type, analysis_type, distance)
  }
  
  # Combine results
  combined_results <- bind_rows(offtarget_results, proximal_results)
  
  if (nrow(combined_results) > 0) {
    # Save individual file results
    output_file <- file.path("offtarget_results", paste0(target_name, "_annotated.csv"))
    write_csv(combined_results, output_file)
    cat("Saved:", output_file, "\n")
  }
  
  return(combined_results)
}

# PROCESS ALL OFF-TARGET FILES FOR SURFACE MARKERS
targets <- c("CD55", "CD81")
all_results <- list()

for (target in targets) {
  # Process individual guides
  guide_files <- list.files("offtarget_data", 
                           pattern = paste0("^", target, "_CRISPRoff_guide[0-9]+\\.csv$"), 
                           full.names = TRUE)
  
  target_results <- list()
  for (i in seq_along(guide_files)) {
    guide_name <- paste0(target, "_CRISPRoff_guide", i)
    result <- process_offtarget_file(guide_files[i], guide_name, target)
    if (!is.null(result)) {
      target_results[[guide_name]] <- result
    }
  }
  
  # Combine all guides for pooled analysis
  if (length(target_results) > 0) {
    # When combining multiple guides, prioritize off-target classification over proximal
    pooled_results <- bind_rows(target_results) %>%
      group_by(gene_symbol) %>%
      summarise(
        target_file = paste0(target, "_pool"),
        locus = dplyr::first(locus),
        entrez_id = dplyr::first(entrez_id),
        mismatch = dplyr::first(mismatch),
        gene_type = dplyr::first(gene_type),
        # If a gene is off-target in ANY guide, classify it as off-target
        analysis_type = if_else(any(analysis_type == "offtarget"), "offtarget", dplyr::first(analysis_type)),
        distance = dplyr::first(distance),
        .groups = "drop"
      )
    
    all_results[[paste0(target, "_pool")]] <- pooled_results
    
    # Save pooled results
    output_file <- file.path("offtarget_results", paste0(target, "_pool_annotated.csv"))
    write_csv(pooled_results, output_file)
    cat("\nSaved pooled results:", output_file, "\n")
  }
}

# STEP 2: CROSS-REFERENCE WITH RNA-SEQ DATA AND CREATE HEATMAPS
cat("\n\nCROSS-REFERENCING WITH RNA-SEQ DATA\n")

# Load RNA-seq results
rna_results <- list()
for (target in targets) {
  file_path <- file.path("results_dge", paste0("DE_", target, "_vs_NTC.csv"))
  if (file.exists(file_path)) {
    rna_results[[target]] <- read_csv(file_path, show_col_types = FALSE)
    cat("Loaded RNA-seq results for", target, "\n")
  }
}

# Create heatmap data for each target pool
for (target in targets) {
  pool_name <- paste0(target, "_pool")
  
  if (!is.null(all_results[[pool_name]]) && !is.null(rna_results[[target]])) {
    cat("\nCreating heatmap data for", target, "...\n")
    
    # Get unique genes from off-target analysis
    analysis_genes <- all_results[[pool_name]] %>%
      dplyr::select(gene_symbol, analysis_type) %>%
      distinct()
    
    # Merge with RNA-seq data
    heatmap_data <- rna_results[[target]] %>%
      inner_join(analysis_genes, by = "gene_symbol") %>%
      mutate(
        significance = adj.P.Val < 0.05
      ) %>%
      # Sort by analysis type with custom order, then by logFC
      mutate(analysis_order = case_when(
        analysis_type == "offtarget" ~ 1,
        analysis_type == "proximal" ~ 2,
        analysis_type == "both" ~ 3,
        analysis_type == "target" ~ 4,
        TRUE ~ 5
      )) %>%
      arrange(analysis_order, desc(abs(logFC))) %>%
      dplyr::select(-analysis_order)
    
    # Ensure target gene is included
    if (!(target %in% heatmap_data$gene_symbol)) {
      target_data <- rna_results[[target]] %>%
        filter(gene_symbol == target) %>%
        mutate(
          analysis_type = "target",
          significance = adj.P.Val < 0.05
        )
      
      if (nrow(target_data) > 0) {
        heatmap_data <- bind_rows(heatmap_data, target_data)
      }
    } else {
      # Update analysis type for target gene
      heatmap_data <- heatmap_data %>%
        mutate(analysis_type = ifelse(gene_symbol == target, "target", analysis_type))
    }
    
    # Save heatmap data
    output_file <- file.path("heatmap_results", paste0(pool_name, "_logFC_data.csv"))
    write_csv(heatmap_data, output_file)
    cat("Saved heatmap data:", output_file, "\n")
    
    # Create guide info
    guide_info <- tibble(
      guide = pool_name,
      target = target,
      num_guides = 3,
      num_offtargets = sum(heatmap_data$analysis_type == "offtarget"),
      num_proximal = sum(heatmap_data$analysis_type == "proximal"),
      num_significant = sum(heatmap_data$significance)
    )
    
    guide_info_file <- file.path("heatmap_results", paste0(pool_name, "_guide_info.csv"))
    write_csv(guide_info, guide_info_file)
    
    # Skip basic heatmap generation - will use enhanced heatmap script instead
  }
}

# Print summary
cat("\n\nANALYSIS SUMMARY\n")
cat("================\n")
for (pool_name in names(all_results)) {
  result_summary <- all_results[[pool_name]] %>%
    group_by(analysis_type) %>%
    summarise(count = n(), .groups = "drop")
  
  cat("\n", pool_name, ":\n", sep = "")
  print(result_summary)
}

cat("\nOFF-TARGET ANALYSIS COMPLETE\n")