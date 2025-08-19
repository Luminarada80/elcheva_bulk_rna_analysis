#! /bin/bash

set -euo pipefail

source activate bioinformatics
module load STAR/2.7.3a
module load fastqc
module load samtools/1.13
module load subread/2.0.6

cd /gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis

SAMPLE_NAME="K562_shCNTL_Rep1"

DATA_DIR="/gpfs/Labs/Uzun/DATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"
FASTQ_DATA_DIR="${DATA_DIR}/FASTQ_FILES/RAW/BATCH_2025_01"

TEST_FILE_R1="${FASTQ_DATA_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
TEST_FILE_R2="${FASTQ_DATA_DIR}/${SAMPLE_NAME}_R2.fastq.gz"

RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"

# ----- Genome Annotation Files -----
genome_index_dir=/gpfs/Labs/Uzun/DATA/GENOMES/INDEX/HUMAN/HG38/STAR
gene_annot_gtf=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf
gene_annot_bed=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf.bed

# --- Trimmomatic paths ---
BASE=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis
TRIMMO_JAR="$BASE/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar"
ADAPTERS="$BASE/trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# ----- Other Dataset Paths -----
METADATA_DIR="/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"
COMPARISON_FILE="${METADATA_DIR}/Comparisons.2025_01.txt"
SAMPLE_ANNOT_FILE="${METADATA_DIR}/Sample_Annotation.2025_01.txt"

SECONDS=0

# Run FastQC
# echo "Running FastQC read quality assessment, saving to ${RESULTS_DIR}/initial_fastqc_reports"
# mkdir -p "${RESULTS_DIR}/initial_fastqc_reports/${SAMPLE_NAME}"
# fastqc $TEST_FILE_R1 -o "${RESULTS_DIR}/initial_fastqc_reports/${SAMPLE_NAME}"

# Adapter Trimming using trimmomatic
# echo "Trimming adapter sequences"
# mkdir -p ${RESULTS_DIR}/trimmed_reads
# java -jar trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -threads 8 -phred33 \
#     $TEST_FILE_R1 \
#     $TEST_FILE_R2 \
#     $RESULTS_DIR/trimmed_reads/${SAMPLE_NAME}_R1_trimmed.fastq \
#     $RESULTS_DIR/trimmed_reads/${SAMPLE_NAME}_R2_trimmed.fastq \
#     ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
#     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# echo "Trimmomatic finished running!"

# Run FastQC on the trimmed reads
# echo "Running FastQC on the trimmed reads"
# mkdir -p "${RESULTS_DIR}/trimming_fastqc_reports"
# fastqc $RESULTS_DIR/trimmed_reads/${SAMPLE_NAME}_R1_trimmed.fastq -o $RESULTS_DIR/trimming_fastqc_reports

# Mapping reads to the genome using STAR
# echo "Running STAR for alignment"
# mkdir -p $RESULTS_DIR/STAR_alignment/${SAMPLE_NAME}/
# STAR --genomeDir "$genome_index_dir" \
#      --runThreadN 8 \
#      --readFilesIn "$RESULTS_DIR/trimmed_reads/${SAMPLE_NAME}_R1_trimmed.fastq" "$RESULTS_DIR/trimmed_reads/${SAMPLE_NAME}_R2_trimmed.fastq" \
#      --readFilesCommand zcat \
#      --twopassMode Basic \
#      --outFileNamePrefix "$RESULTS_DIR/STAR_alignment/${SAMPLE_NAME}/${SAMPLE_NAME}_" \
#      --outSAMtype BAM SortedByCoordinate \
#      --limitBAMsortRAM 16000000000 \
#      --outSJfilterReads Unique \
#      --outSAMattributes Standard \
#      --outSAMattrRGline ID:$SAMPLE_NAME SM:$SAMPLE_NAME PL:ILLUMINA

# echo "STAR finished running!"

# Build a BAM index
echo "Building BAM index"
module load samtools/1.13
samtools index ${RESULTS_DIR}/STAR_alignment/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam

# Filter unique reads
echo "Filtering unique reads"
mkdir -p ${RESULTS_DIR}/TEMP_SORTING/test/
samtools view -bh -q 255 ${RESULTS_DIR}/STAR_alignment/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam \
    > ${RESULTS_DIR}/STAR_alignment/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.bam

# Compute FPKM with Cufflinks
echo "Running Cufflinks for FPKM"
module load cufflinks/2.2.1
mkdir -p ${RESULTS_DIR}/cufflinks_fpkm_results/
cufflinks-2.2.1.Linux_x86_64/cufflinks -p 8 \
    --library-type fr-firststrand \
    -o ${RESULTS_DIR}/cufflinks_fpkm_results/ \
    -G $gene_annot_gtf \
    ${RESULTS_DIR}/STAR_alignment/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.bam

echo "FPKM computed"

# Compute integer read counts with featureCounts
echo "Running featureCounts for read counts"
mkdir -p "${RESULTS_DIR}/featurecounts/${SAMPLE_NAME}"
featureCounts \
    -T 8 \
    -p \
    -B \
    -C \
    -s 2 \
    -a "$gene_annot_gtf" \
    -t exon \
    -g gene_id \
    -o "${RESULTS_DIR}/featurecounts/${SAMPLE_NAME}/${SAMPLE_NAME}_gene_counts.txt" \
    "${RESULTS_DIR}/STAR_alignment/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.bam"

echo "featureCounts done."

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."