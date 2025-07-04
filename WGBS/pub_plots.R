library(tidyverse)
library(ggplot2)
library(GenomicRanges)

# Define custom colors
positive_color <- "#E41A1C"  # A beautiful color for positive values
negative_color <- "#377EB8"  # A beautiful color for negative values

gtf <- rtracklayer::import('annotations/hg38.ncbiRefSeq.gtf.gz')
gene_annotation=as.data.frame(gtf)

generate_plot <- function(file_name, sample_name) {
  DMRs_ann <- read_csv(file_name) %>%
    mutate(CHR=chr,BP=end+start/2) %>% 
    mutate(P=sign(meth.diff) * -log10(pvalue))
  
  DMRs_ann$CHR <- gsub("chr(\\d+)", "chr\\1", DMRs_ann$CHR)
  DMRs_ann$CHR <- gsub("chr([0-9])$", "chr0\\1", DMRs_ann$CHR)

  don <- DMRs_ann %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(DMRs_ann, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    mutate( is_annotate=ifelse( P>40, "yes", "no")) 
  
  # #Convert your data and the gene annotation data to GRanges objects
  # gr_mehtylation_data <- GRanges(seqnames = DMRs_ann$CHR, ranges = IRanges(start = DMRs_ann$start, end = DMRs_ann$end))
  # gr_gene_annotation <- GRanges(seqnames = gene_annotation$seqnames, ranges = IRanges(start = gene_annotation$start, end = gene_annotation$end), gene = gene_annotation$gene_name)
  # 
  # # Find overlaps
  # overlap_counts <- countOverlaps(gr_mehtylation_data, gr_gene_annotation)
  # 
  # # Find indices of ranges with multiple overlaps
  # multiple_overlap_indices <- which(overlap_counts > 1)
  # 
  # # For each range with multiple overlaps, assign the gene with the highest overlap
  # for (i in multiple_overlap_indices) {
  #   overlaps_for_range <- findOverlaps(gr_mehtylation_data[i], gr_gene_annotation)
  #   overlap_lengths <- width(intersect(gr_mehtylation_data[i], gr_gene_annotation[subjectHits(overlaps_for_range)]))
  #   gene_with_max_overlap <- as.character(gr_gene_annotation$gene[subjectHits(overlaps_for_range)][which.max(overlap_lengths)])
  #   DMRs_ann$gene[i] <- gene_with_max_overlap
  # }

  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2, max=max(BPcum))

  don$abs_P <- abs(don$P)  # Calculate the absolute value of P

  manhattan_plot <- ggplot(don, aes(x=BPcum, y=P)) +
    
    # Show all points
    geom_point(aes(color = ifelse(P > 0, "Positive", "Negative"), size = abs_P/10), 
               # alpha = ifelse(don$abs_P > 2, .1, 0.01)
               alpha = ifelse(don$abs_P <= 2, 0.01, (don$abs_P - 1) / max(don$abs_P - 1))
    ) +
    scale_color_manual(values = c("Positive" = positive_color, "Negative" = negative_color)) +
    scale_size_continuous(range = c(.5, .5)) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks=axisdf$center ) +
    scale_y_continuous( expand = c(0, 0), limits = c(-100, 150)) +     # remove space between plot area and x axis
    
    geom_vline(xintercept = axisdf$max, linetype='dotted', linewidth = .3) +
    
    # Add highlighted points
    # geom_point(data=subset(don, is_highlight=="yes"), color="red", alpha=0.4, size=.81) +
    
    xlab('Genomic Location') +
    ylab('-log10(pvalue) * sign(methylation difference)') +
    ggtitle(paste0(sample_name)) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      plot.title = element_text(face = "bold", size = rel(0.6), hjust = 0.5), 
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA), 
      # panel.border = element_rect(colour = NA), 
      axis.title = element_text(face = "bold", size = rel(1)), 
      axis.title.y = element_text(angle = 90, vjust = 2), 
      axis.title.x = element_text(vjust = -0.2), 
      axis.text.x = element_text(vjust = -0.1, angle = 90),
      
      plot.margin = unit(c(10, 5, 5, 5), "mm"), 
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"), 
      strip.text = element_text(face = "bold"),
      
      legend.position="none",
      # panel.border = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      # axis.line = element_line(linewidth = .1)
      
    )

  ggsave(paste0('plots/', sample_name, '.png'), manhattan_plot, width = 30, height=7.5, units = 'cm')
}

generate_plot('results/CD55_methylkit_1000_100.csv', 'CD55')
generate_plot('results/CD81_methylkit_1000_100.csv', 'CD81')
generate_plot('results/EMPTY_methylkit_1000_100.csv', 'EMPTY Electroporation')