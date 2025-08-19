#!/bin/bash -l

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail

source activate bioinformatics
module load STAR/2.7.3a
module load fastqc
module load cufflinks/2.2.1
module load samtools/1.13
module load subread/2.0.6

BASE=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis
cd $BASE

DATA_DIR="/gpfs/Labs/Uzun/DATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"
FASTQ_DATA_DIR="${DATA_DIR}/FASTQ_FILES/RAW/BATCH_2025_01"

TEST_FILE_R1="${FASTQ_DATA_DIR}/${SAMPLE_NAME}_R1.fastq.gz" # Only use R1

RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/${SAMPLE_NAME}" #
OUTPUT_DIR="${RESULTS_DIR}/read_alignment"
mkdir -p $OUTPUT_DIR

# ----- Genome Annotation Files -----
genome_index_dir=/gpfs/Labs/Uzun/DATA/GENOMES/INDEX/HUMAN/HG38/STAR
gene_annot_gtf=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf
gene_annot_bed=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf.bed


# ----- Other Dataset Paths -----
METADATA_DIR="/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"
COMPARISON_FILE="${METADATA_DIR}/Comparisons.2025_01.txt"
SAMPLE_ANNOT_FILE="${METADATA_DIR}/Sample_Annotation.2025_01.txt"

TRIMMO_JAR="$BASE/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar"
ADAPTERS="$BASE/trimmomatic-0.39/adapters/TruSeq3-SE.fa"

SECONDS=0


# This pipeline runs QC and alignment for reads generated using Lexogen's QuantSeq 3' mRNA-Seq
# V2 Library Prep Kit Forward with 12 nucleotide Unique Dual Indices (UDIs).

# Notes:
#  1. Only R1 reads are used, as R2 reads start with poly(T) and are low quality

echo "Output Parent Directory: ${OUTPUT_DIR}"

# Run FastQC
echo "Running FastQC read quality assessment"
echo "    - Saving to Output Directory/01_raw_read_fastqc"
mkdir -p "${OUTPUT_DIR}/01_raw_read_fastqc"
fastqc $TEST_FILE_R1 -o "${OUTPUT_DIR}/01_raw_read_fastqc"

# Adapter trimming using trimmomatic
echo "Trimming adapters and low-quality reads using Trimmomatic"
echo "    - Saving to Output Directory/02_trimmomatic"
mkdir -p "${OUTPUT_DIR}/02_trimmomatic"
java -jar $TRIMMO_JAR SE -threads 8 -phred33 \
    $TEST_FILE_R1 ${OUTPUT_DIR}/02_trimmomatic/${SAMPLE_NAME}_R1_trimmed.fastq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

# Run FastQC on the trimmed reads
echo "Running FastQC on the trimmed reads"
echo "    - Saving to Output Directory/03_trimmed_read_fastqc"
mkdir -p "${OUTPUT_DIR}/03_trimmed_read_fastqc"
fastqc ${OUTPUT_DIR}/02_trimmomatic/${SAMPLE_NAME}_R1_trimmed.fastq.gz -o ${OUTPUT_DIR}/03_trimmed_read_fastqc

# Mapping reads to the genome using STAR
echo "Running STAR for alignment"
echo "    - Saving to Output Directory/04_star_alignment"
mkdir -p ${OUTPUT_DIR}/04_star_alignment
STAR --genomeDir "$genome_index_dir" \
     --runThreadN 8 \
     --readFilesIn "${OUTPUT_DIR}/02_trimmomatic/${SAMPLE_NAME}_R1_trimmed.fastq.gz" \
     --readFilesCommand zcat \
     --twopassMode Basic \
     --outFileNamePrefix "${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_" \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 16000000000 \
     --outSJfilterReads Unique \
     --outSAMattributes Standard \
     --outSAMattrRGline ID:$SAMPLE_NAME SM:$SAMPLE_NAME PL:ILLUMINA

echo "      STAR finished running!"

# Build a BAM index
echo "Building BAM index"
samtools index ${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam

# Filter unique reads
echo "Filtering unique reads"
echo "    - Saving to Output Directory/04_star_alignment/${SAMPLE_NAME}_Aligned.bam"
samtools view -bh -q 255 ${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam \
    > ${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_Aligned.bam

# Compute FPKM with Cufflinks
echo "Running Cufflinks for FPKM"
echo "    - Saving to Output Directory/05_cufflinks_fpkm"

mkdir -p ${OUTPUT_DIR}/05_cufflinks_fpkm
cufflinks -p 8 \
    --library-type fr-secondstrand \
    -o ${OUTPUT_DIR}/05_cufflinks_fpkm/ \
    -G $gene_annot_gtf \
    ${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_Aligned.bam

echo "      FPKM computed"

# Compute integer read counts with featureCounts
echo "Running featureCounts for read counts"
echo "    - Saving to Output Directory/06_featureCounts"

mkdir -p "${OUTPUT_DIR}/06_featureCounts"
featureCounts \
    -T 8 \
    -s 1 \
    -a "$gene_annot_gtf" \
    -t exon \
    -g gene_id \
    -o "${OUTPUT_DIR}/06_featureCounts/${SAMPLE_NAME}_gene_counts.txt" \
    "${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_Aligned.bam"

mkdir -p ${BASE}/gene_counts
mv "${OUTPUT_DIR}/06_featureCounts/${SAMPLE_NAME}_gene_counts.txt" "${RESULTS_DIR}/gene_counts/${SAMPLE_NAME}_gene_counts.txt"

echo "featureCounts done."

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."