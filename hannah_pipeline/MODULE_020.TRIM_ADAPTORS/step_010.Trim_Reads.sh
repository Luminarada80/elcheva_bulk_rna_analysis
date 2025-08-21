#!/bin/bash -l
#SBATCH --job-name=Trim
#SBATCH --output=Trim_results.%j.%N.txt
#SBATCH --error=Trim_errors.%j.%N.err
#SBATCH -p compute
#SBATCH -A lab_uzun
#SBATCH --nodes 2
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=50G
srun esm5360

project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

raw_fastq_dir=/gpfs/Labs/Uzun/DATA/PROJECTS/$project_name/FASTQ_FILES/RAW/BATCH_$batch_name/
trim_fastq_dir=/gpfs/Labs/Uzun/DATA/PROJECTS/$project_name/FASTQ_FILES.TRIMMED/$batch_name/

mkdir $trim_fastq_dir

threads=4
adapters=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/$project_name/MODULE_020.TRIM_ADAPTERS/CURRENT/TruSeq3_PE.fa

cd $raw_fastq_dir

module load trimmomatic/0.38

for sample in ${raw_fastq_dir}/*_R1.fastq.gz;
do

base_name=$(basename $sample _R1.fastq.gz)


forward_input=${raw_fastq_dir}/${base_name}_R1.fastq.gz
reverse_input=${raw_fastq_dir}/${base_name}_R2.fastq.gz


forward_paired_output=${trim_fastq_dir}/${base_name}_R1_paired.fastq.gz
forward_unpaired_output=${trim_fastq_dir}/${base_name}_R1_unpaired.fastq.gz
reverse_paired_output=${trim_fastq_dir}/${base_name}_R2_paired.fastq.gz
reverse_unpaired_output=${trim_fastq_dir}/${base_name}_R2_unpaired.fastq.gz


 trimmomatic PE -threads $threads -phred33 \
        $forward_input $reverse_input \
        $forward_paired_output $forward_unpaired_output \
        $reverse_paired_output $reverse_unpaired_output \
        ILLUMINACLIP:${adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
echo "Trimmomatic processing complete."
