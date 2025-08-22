#!/bin/bash -l
#SBATCH --job-name=Trim
#SBATCH --output=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/LOGS/trimming/Trim_results.%j.log
#SBATCH --error=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/LOGS/trimming/Trim_errors.%j.err
#SBATCH -p compute
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --mem=64G

project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

raw_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/RAW/BATCH_${batch_name}
trim_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/FASTQ_FILES.TRIMMED/BATCH_${batch_name}

BASE=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/${project_name}/elcheva_bulk_rna_analysis

mkdir -p $trim_fastq_dir

threads=4

cd $raw_fastq_dir

# module load trimmomatic/0.38

for sample in ${raw_fastq_dir}/*_R1.fastq.gz;
do

base_name=$(basename "$sample" _R1.fastq.gz)


forward_input=${raw_fastq_dir}/${base_name}_R1.fastq.gz
reverse_input=${raw_fastq_dir}/${base_name}_R2.fastq.gz


forward_paired_output=${trim_fastq_dir}/${base_name}_R1_paired.fastq.gz
forward_unpaired_output=${trim_fastq_dir}/${base_name}_R1_unpaired.fastq.gz
reverse_paired_output=${trim_fastq_dir}/${base_name}_R2_paired.fastq.gz
reverse_unpaired_output=${trim_fastq_dir}/${base_name}_R2_unpaired.fastq.gz

TRIMMO_JAR="$BASE/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar"
ADAPTERS="$BASE/trimmomatic-0.39/adapters/TruSeq3-SE.fa"


java -jar $TRIMMO_JAR PE -threads $threads -phred33 \
        $forward_input $reverse_input \
        $forward_paired_output $forward_unpaired_output \
        $reverse_paired_output $reverse_unpaired_output \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
echo "Trimmomatic processing complete."
