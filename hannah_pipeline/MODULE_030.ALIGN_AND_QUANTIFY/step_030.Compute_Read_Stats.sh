project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

meta_data_dir=/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/
batch_sample_list_file=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/sample_names_list.txt
batch_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/FASTQ_FILES.TRIMMED/BATCH_${batch_name}
batch_bam_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/BAM_FILES/$batch_name/


batch_read_stats_dir=/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/$project_name/READ_STATS/${batch_name}/INDIVIDUAL_SAMPLES/

genome_index_dir=/gpfs/Labs/Uzun/DATA/GENOMES/INDEX/HUMAN/HG38/STAR/
gene_annot_gtf=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf
gene_annot_bed=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf.bed

batch_job_dir=/gpfs/Labs/Uzun/JOBS/PROJECTS/${project_name}/ALIGNMENT/

mkdir -p $batch_job_dir

mkdir -p $batch_bam_dir
mkdir -p $batch_bam_dir/TEMP_SORTING # for sorting
mkdir -p $batch_fpkm_dir
mkdir -p $batch_read_count_dir
mkdir -p $batch_read_stats_dir

IFSs=$'\n' read -d '' -r -a sample_names < $batch_sample_list_file
len=${#sample_names[@]}

echo "There are $len samples"
echo "The sample names are:"
echo ${sample_names[@]}


echo $batch_bam_dir/
echo $batch_job_dir
ls $batch_job_dir

for sample_name in ${sample_names[@]}
do
  echo $sample_name
  
  job_file=$batch_job_dir/read_stats.${sample_name}.sh
  log_file=$batch_job_dir/read_stats.${sample_name}.log
  err_file=$batch_job_dir/read_stats.${sample_name}.log
  rm -f $log_file
  rm -f $err_file

  library=$sample_name
  echo $job_file
  echo $log_file
  echo "#! /bin/bash" > $job_file
  echo "#SBATCH --output=$log_file">> $job_file
  echo "#SBATCH --error=$err_file">> $job_file
  echo "#SBATCH --ntasks=1">> $job_file
  echo "#SBATCH --time=2-0">> $job_file
  echo "#SBATCH --partition=compute">> $job_file
  echo "#SBATCH --account=lab_uzun">> $job_file
  echo "#SBATCH --cpus-per-task=9">> $job_file
  echo "#SBATCH --mem-per-cpu=20G">> $job_file
  #echo "#SBATCH --mem-per-cpu=10G">> $job_file
  echo ".  /etc/profile.d/modules.sh">> $job_file

  echo "source activate bioinformatics">> $job_file

  # Filter unique reads
  echo "cd $batch_bam_dir">> $job_file
  echo "module load samtools/1.13">> $job_file
  echo "mkdir -p TEMP_SORTING/${sample_name}/">> $job_file	
  echo "samtools sort -T TEMP_SORTING/${sample_name}/ -o ${sample_name}_Aligned.sorted.bam ${sample_name}_Aligned.bam">> $job_file
  echo "samtools index ${sample_name}_Aligned.sorted.bam">> $job_file
  echo "echo Reads filtered.">> $job_file

  # Compute read (intron/exon/intergenic) statistics 
  echo "echo Computing read stats.">> $job_file
  echo "cd $batch_bam_dir">> $job_file
  # echo "module load python/2.7.15">> $job_file
  echo "python /gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/MODULE_030.ALIGN_AND_QUANTIFY/ReadStats.py $batch_bam_dir/${sample_name}_Aligned.sorted.bam $gene_annot_bed $batch_read_stats_dir/${sample_name}.Read_stats.txt  ">> $job_file
  echo "echo Read stats computed.">> $job_file
  sbatch $job_file

done

sleep 5
squeue -u hqv5055


echo $batch_bam_dir
ls $batch_bam_dir/
echo $batch_job_dir
ls $batch_job_dir* | head





