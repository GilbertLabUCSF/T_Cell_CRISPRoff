fastqDIR='bam_files/deduplicated/'

# for bed in ${fastqDIR}/*.bismark.cov.gz; do
for bam in ${fastqDIR}/*bismark_bt2_pe.srt.bam; do
    out=`basename $bam`; out=${out/_R1_val_1_bismark_bt2_pe/_qc_hg38};out=${out/.gz/}
    echo $out
    bedtools intersect -wa -b qc_genes_hg38.bed -a ${bam} > igv/exp3/$out
    echo "Done!"    
done