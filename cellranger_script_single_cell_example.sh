/volumes/seq/code/3rd_party/cellranger/cellranger-3.1.0/cellranger count --id sample_id \
    --fastqs /volumes/seq/projects/HBCA/HBCA_10x_RNA/fastq_output/sample_id  \
    --sample  sample_id \
    --transcriptome /volumes/seq/code/3rd_party/10X/refdata-cellranger-GRCh38-3.0.0 
    --localcores=20 \
    --localmem=300