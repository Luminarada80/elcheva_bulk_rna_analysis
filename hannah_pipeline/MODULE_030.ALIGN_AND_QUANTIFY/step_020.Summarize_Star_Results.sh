project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

meta_data_dir=/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/
batch_sample_list_file=/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis/hannah_pipeline/sample_names_list.txt
batch_fastq_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/${project_name}/FASTQ_FILES/FASTQ_FILES.TRIMMED/BATCH_${batch_name}
batch_bam_dir=/gpfs/Labs/Uzun/DATA/GRANT_APPS/$project_name/BULK_RNA_SEQ/BAM_FILES/$batch_name/

result_dir=/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/$project_name/STAR_ALIGNMENT_STATISTICS/
mkdir -p $result_dir
####


cd $batch_bam_dir 
find ./ -name '*Log.final.out*' | sort -n | xargs awk 'FNR==1 {print "XXXXXXXXXXXXXXXXXX " FILENAME " XXXXXXXXXXXXXXXXXX"}{print}' > $result_dir/Summary.STAR.txt

cd $result_dir 

stat_file_name=STAR_Alignment_Statistics.${batch_name}.txt

grep XXXXXXXXXXXXXXXXXX Summary.STAR.txt > temp1
grep "Number of input reads" Summary.STAR.txt | awk '{print $NF}' > temp2
grep "Uniquely mapped reads number" Summary.STAR.txt | awk '{print $NF}'  > temp3
grep "Uniquely mapped reads %" Summary.STAR.txt | awk '{print $NF}'  > temp4

echo "Sample Library_Size Uniquely_Aligned_Read_Count Unique_Alignment_Rate" | sed s/" "/"\t"/g > $stat_file_name
paste temp1 temp2 temp3 temp4 | sed s/"XXXXXXXXXXXXXXXXXX ."/""/g  | sed s/"_Log.final.out XXXXXXXXXXXXXXXXXX"/""/g | sed s/"\/batch_2_"/""/  | sed s/"\/"/""/ >> $stat_file_name

cat $stat_file_name







