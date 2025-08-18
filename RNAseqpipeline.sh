#! /bin/bash

source activate bioinformatics

module load STAR/2.7.3a
module load fastqc

cd /gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis

DATA_DIR="/gpfs/Labs/Uzun/DATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"
METADATA_DIR="/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"

FASTQ_DATA_DIR="${DATA_DIR}/FASTQ_FILES/RAW/BATCH_2025_01"

COMPARISON_FILE="${METADATA_DIR}/Comparisons.2025_01.txt"
SAMPLE_ANNOT_FILE="${METADATA_DIR}/Sample_Annotation.2025_01.txt"

GENOME_FASTA=/gpfs/Labs/Uzun/DATA/GENOMES/SEQUENCE/HUMAN/HG38/hg38.fa
genome_index_dir=/gpfs/Labs/Uzun/DATA/GENOMES/INDEX/HUMAN/HG38/STAR
gene_annot_gtf=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf
gene_annot_bed=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf.bed

RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/"

TEST_FILE_R1="${FASTQ_DATA_DIR}/K562_shCNTL_Rep1_R1.fastq.gz"
TEST_FILE_R2="${FASTQ_DATA_DIR}/K562_shCNTL_Rep1_R2.fastq.gz"

SECONDS=0

# STEP 1: Run FastQC
# echo "Running FastQC read quality assessment, saving to ${RESULTS_DIR}/initial_fastqc_reports"
# fastqc $TEST_FILE -o "${RESULTS_DIR}/initial_fastqc_reports"

# Adapter Trimming using trimmomatic
# echo "Trimming adapter sequences"
# java -jar trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -threads 8 -phred33 \
#     $TEST_FILE_R1 \
#     $TEST_FILE_R2 \
#     $RESULTS_DIR/K562_shCNTL_Rep1_R1_trimmed.fastq \
#     $RESULTS_DIR/K562_shCNTL_Rep1_R2_trimmed.fastq \
#     ILLUMINACLIP:"$ADAPTERS":2:30:10 \
#     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# echo "Trimmomatic finished running!"

# Trimming QC
# fastqc $RESULTS_DIR/K562_shCNTL_Rep1_R1_trimmed.fastq -o $RESULTS_DIR/trimming_fastqc_reports

# # Mapping reads to the genome using STAR
# echo "Running STAR for alignment"
# STAR --genomeDir "$genome_index_dir" \
#      --runThreadN 8 \
#      --readFilesIn "$TEST_FILE_R1" "$TEST_FILE_R2" \
#      --readFilesCommand zcat \
#      --twopassMode Basic \
#      --outFileNamePrefix "$RESULTS_DIR/K562_shCNTL_" \
#      --outSAMtype BAM SortedByCoordinate \
#      --limitBAMsortRAM 16000000000 \
#      --outSJfilterReads Unique \
#      --outSAMattributes Standard \
#      --outSAMattrRGline ID:K562_shCNTL SM:K562_shCNTL PL:ILLUMINA

# echo "STAR finished running!"

# # Build a BAM index
# echo "Building BAM index"
# module load samtools/1.13
# samtools index ${RESULTS_DIR}/K562_shCNTL_Aligned.sortedByCoord.out.bam

# # Filter unique reads
# echo "Filtering unique reads"
# mkdir -p ${RESULTS_DIR}/TEMP_SORTING/test/
# samtools view -bh -q 255 ${RESULTS_DIR}/K562_shCNTL_Aligned.sortedByCoord.out.bam \
#     > ${RESULTS_DIR}/K562_shCNTL_Aligned.bam

# Compute FPKM with Cufflinks
echo "Running Cufflinks for FPKM"
module load cufflinks/2.2.1
cufflinks-2.2.1.Linux_x86_64/cufflinks -p 8 \
    --library-type fr-firststrand \
    -o $RESULTS_DIR/cufflinks_fpkm_results/ \
    -G $gene_annot_gtf \
    $RESULTS_DIR/K562_shCNTL_Aligned.bam

echo "FPKM computed"

# Compute integer read counts with featureCounts
echo "Running featureCounts for read counts"


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."