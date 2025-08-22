project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

batch_fpkm_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/FPKM/${batch_name}/INDIVIDUAL_SAMPLES/
batch_read_count_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/READ_COUNTS/${batch_name}/INDIVIDUAL_SAMPLES/


batch_matrix_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/${batch_name}/
mkdir -p $batch_matrix_dir

module load R

# Rscript  /gpfs/Labs/Uzun/SCRIPTS/COMMON/RNA_SEQ/CURRENT/merge_fpkm.R  $batch_fpkm_dir $batch_matrix_dir
# ls -lh $batch_matrix_dir

Rscript /gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/MODULE_030.ALIGN_AND_QUANTIFY/merge_read_counts.R $batch_read_count_dir/GENE_LEVEL  $batch_matrix_dir/Read_Count_Matrix.All_genes
ls -lh $batch_matrix_dir

Rscript /gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/MODULE_030.ALIGN_AND_QUANTIFY/merge_read_counts.R $batch_read_count_dir/TRANSCRIPT_LEVEL  $batch_matrix_dir/Read_Count_Matrix.All_transcripts
ls -lh $batch_matrix_dir
