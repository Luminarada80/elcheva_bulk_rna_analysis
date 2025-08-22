project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

meta_data_dir=/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/
batch_sample_list_file=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/sample_names_list.txt
batch_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/FASTQ_FILES.TRIMMED/BATCH_${batch_name}
batch_bam_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/BAM_FILES/$batch_name/

batch_fpkm_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/FPKM/${batch_name}/INDIVIDUAL_SAMPLES/
batch_read_count_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/READ_COUNTS/${batch_name}/INDIVIDUAL_SAMPLES/

genome_index_dir=/gpfs/Labs/Uzun/DATA/GENOMES/INDEX/HUMAN/HG38/STAR/
gene_annot_gtf=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf
gene_annot_bed=/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf.bed

batch_job_dir=/gpfs/Labs/Uzun/JOBS/PROJECTS/${project_name}/ALIGNMENT/

mkdir -p $batch_job_dir

mkdir -p $batch_bam_dir
mkdir -p $batch_bam_dir/TEMP_SORTING # for sorting
mkdir -p $batch_fpkm_dir
mkdir -p $batch_read_count_dir
# mkdir -p $batch_read_stats_dir

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
  
  job_file=$batch_job_dir/align.${sample_name}.sh
  log_file=$batch_job_dir/align.${sample_name}.log
  err_file=$batch_job_dir/align.${sample_name}.log
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
  # echo "srun hostname">> $job_file
  # echo ".  /etc/profile.d/modules.sh">> $job_file

  # #STAR Align
  # echo "echo Running STAR for alignment.">> $job_file
  # echo "cd $batch_fastq_dir">> $job_file
  # echo "module load STAR/2.7.3a">> $job_file
  # echo "STAR --genomeDir $genome_index_dir \\">> $job_file
  # echo "--runThreadN 8 \\">> $job_file
  # echo "--readFilesIn ${sample_name}_R1_paired.fastq.gz ${sample_name}_R2_paired.fastq.gz \\">> $job_file
  # echo "--outFileNamePrefix $batch_bam_dir/${sample_name}_ \\">> $job_file
  # echo "--outSAMtype BAM SortedByCoordinate \\">> $job_file
  # echo "--outSJfilterReads Unique \\">> $job_file
  # echo "--readFilesCommand zcat \\">> $job_file
  # echo "--outSAMattributes Standard ">> $job_file
  # echo "echo Alignment completed.">> $job_file

  # #BAM index
  # echo "echo Building BAM index.">> $job_file
  # echo "cd $batch_bam_dir">> $job_file
  # echo "module load samtools/1.10.0">> $job_file
  # echo "samtools index ${sample_name}_Aligned.sortedByCoord.out.bam">> $job_file
  # echo "echo Index built.">> $job_file

  # # Filter unique reads
  # echo "cd $batch_bam_dir">> $job_file
  # echo "mkdir -p TEMP_SORTING/${sample_name}/">> $job_file
  # echo "module load samtools/1.10.0">> $job_file
  # echo "samtools view -bh -q 255 ${sample_name}_Aligned.sortedByCoord.out.bam > ${sample_name}_Aligned.bam" >> $job_file	

  # #Compute FPKM with Cufflinks
  # mkdir -p $batch_fpkm_dir/$sample_name
  # echo "echo Running Cufflinks for FPKM.">> $job_file
  # echo "module load cufflinks/2.2.1">> $job_file
  # echo "rm -r $batch_fpkm_dir/$sample_name/">> $job_file
  # echo "cufflinks -p 8 --library-type fr-firststrand  -o $batch_fpkm_dir/$sample_name/  -G $gene_annot_gtf $batch_bam_dir/${sample_name}_Aligned.bam" >> $job_file
  # echo "python /gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/MODULE_030.ALIGN_AND_QUANTIFY/Select_Greater_FPKM.py $batch_fpkm_dir/$sample_name/genes.fpkm_tracking > $batch_fpkm_dir/$sample_name/genes.fpkm_tracking.max.txt" >> $job_file
  # echo "echo FPKM computed.">> $job_file

  # Compute integer read counts with featureCounts (CLI from Subread)
  echo "echo 'Running featureCounts (CLI) for read counts.'" >> $job_file
  echo "module load subread/2.0.6" >> $job_file   # or conda env with subread installed

  # locate BAM (adjust name if needed, STAR usually outputs *_Aligned.sortedByCoord.out.bam)
  echo "BAM_FILE=$batch_bam_dir/${sample_name}_Aligned.sortedByCoord.out.bam" >> $job_file
  echo "if [[ ! -s \"\$BAM_FILE\" ]]; then echo \"[ERROR] BAM not found: \$BAM_FILE\"; exit 1; fi" >> $job_file

  # gene-level counts
  echo "mkdir -p $batch_read_count_dir/GENE_LEVEL" >> $job_file
  echo "featureCounts -T 8 -p -s 0 -a $gene_annot_gtf -t exon -g gene_id -o $batch_read_count_dir/GENE_LEVEL/${sample_name}.gene.txt \$BAM_FILE" >> $job_file

  # transcript-level counts
  echo "mkdir -p $batch_read_count_dir/TRANSCRIPT_LEVEL" >> $job_file
  echo "featureCounts -T 8 -p -s 0 -a $gene_annot_gtf -t exon -g transcript_id -o $batch_read_count_dir/TRANSCRIPT_LEVEL/${sample_name}.tx.txt \$BAM_FILE" >> $job_file

  echo "echo Read count computed." >> $job_file



  sbatch $job_file

done

sleep 5
squeue -u esm5360


echo $batch_bam_dir
ls $batch_bam_dir/
echo $batch_job_dir
ls $batch_job_dir* | head





