# Co065 Analysis Pipeline

Run the scripts in the following order:

1. `Rscript src/01_Co065_rna_seq_analysis.R`
   - RNA-seq differential expression analysis using limma/edgeR
   - Creates DE results for each target (KO and KD)
   - Generates volcano plots and GO enrichment analysis
   - Outputs to `results_Co065/` and `plots_dge/`

2. `Rscript src/02_Co065_unified_offtarget_analysis.R`
   - Processes off-target data files
   - Creates heatmaps for each guide
   - Generates summary tables
   - Outputs to `Co065_offtarget_results/` and `Co065_heatmap_results/`

3. `Rscript src/03_Co065_enhanced_volcano_plots.R`
   - Creates enhanced volcano plots using DE results from step 1
   - Highlights genes significant in both KD and KO
   - Generates combined plots
   - Outputs to `Co065_enhanced_volcano_plots/`

4. `Rscript src/04_Co065_enhanced_heatmaps.R`
   - Creates enhanced heatmaps with KD∩KO highlighting
   - Purple tracks show genes significant in both conditions
   - Generates summary report
   - Outputs to `Co065_enhanced_heatmaps/`

## Output Directories
- `results_Co065/` - RNA-seq DE results and GO enrichment
- `plots_dge/` - Differential expression plots
- `Co065_offtarget_results/` - Off-target gene annotations
- `Co065_heatmap_results/` - Individual guide heatmaps & data
- `Co065_enhanced_volcano_plots/` - Enhanced volcano plots
- `Co065_enhanced_heatmaps/` - Enhanced heatmaps

## Required Input Files
- `rna_seq_meta_key.csv` - Sample metadata
- `data/tximport_object.rds` - Transcript abundance data
- `data/dge_object.rds` - DGE object
- `data/normalized_counts.csv` - Normalized count matrix
- `offtarget_data/` - Directory with IDT off-target CSV files