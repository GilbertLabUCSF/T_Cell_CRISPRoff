# Integrated Epigenetic and Genetic Programming of Primary Human T Cells

## Overview
This repository contains the code and data analysis for the manuscript "Integrated Epigenetic and Genetic Programming of Primary Human T Cells" by Goudy et al.

## Project Structure

### RNA-seq Analysis
- **Co065_multiple_donors/**: RNA-seq data from multiple donors comparing CRISPRoff vs CRISPR-Cas9 knockout
- **CRISPRoff_KO_MED12/**: Comparison of MED12 silencing between CRISPRoff and CRISPR-Cas9
- **CRISPRoff_surface_markers/**: Comparison of different Cas9 effector proteins on surface marker expression

### Methylation Analysis
- **laine_methylation_recalc/**: Bisulfite sequencing analysis for methylation profiling

## Key Analyses

### 1. CRISPRoff vs CRISPR-Cas9 Comparison
Demonstrates comparable gene silencing efficiency between epigenetic (CRISPRoff) and genetic (CRISPR-Cas9) approaches.

### 2. Multi-target Analysis
Evaluates CRISPRoff performance across multiple therapeutically relevant genes including:
- FAS
- MED12
- PTPN2
- RASA2
- RC3H1
- SUV39H1

### 3. Surface Marker Analysis
Compares different Cas9 variants for targeting T cell surface markers:
- CD55
- CD81
- CD151

## Dependencies
- R (v4.3+)
- DESeq2
- edgeR
- methylKit
- bsseq
- ggplot2
- tidyverse

## Data Availability
Raw sequencing data and processed files are organized in respective subdirectories.

## Contact
For questions regarding this repository, please contact the corresponding authors:
- justin.eyquem@ucsf.edu
- brian.shy@ucsf.edu
- alexander.marson@ucsf.edu
- luke@arcinstitute.org