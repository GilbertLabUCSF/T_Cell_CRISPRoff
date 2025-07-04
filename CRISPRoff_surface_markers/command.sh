nextflow run nf-core/rnaseq \
    --input samplesheet_nfcore.csv \
    --outdir results \
    --gtf /large_storage/gilbertlab/refdata/bulkRNA_seq/GRCh38/gencode.v47.basic.annotation.gtf.gz \
    --fasta /large_storage/gilbertlab/refdata/bulkRNA_seq/GRCh38/GRCh38.p14.genome.fa.gz \
    --igenomes_ignore \
    --genome 'null' \
    -profile apptainer \
    --max_memory '256 GB' \
    -resume