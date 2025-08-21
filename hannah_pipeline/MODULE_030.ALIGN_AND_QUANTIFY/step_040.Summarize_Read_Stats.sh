project_name=2025.NYNRIN.DOD.ELCHEVA
batch_name="2025_01"

individual_read_stats_dir=/gpfs/Labs/Uzun/RESULTS/PROJECTS/$project_name/READ_STATS/${batch_name}/INDIVIDUAL_SAMPLES/
merged_read_stats_dir=/gpfs/Labs/Uzun/RESULTS/PROJECTS/$project_name/READ_STATS/MERGED/

mkdir -p $merged_read_stats_dir

current=$PWD

cd $individual_read_stats_dir

echo "Sample " > results_01.txt
echo "Exon " >> results_01.txt
echo "Intron " >> results_01.txt
echo "Junction " >> results_01.txt
echo "Intergenic " >> results_01.txt
echo "Unique " >> results_01.txt


for filename in $(find -name '*Read_stats*txt' | sort -n)
do

    samplename=$filename
    echo $samplename;
    cat results_01.txt > output_temp
    awk '  NR==1 { print "Sample name : " FILENAME } {if (NR>1) {print}}' $filename > track_temp.txt
    cut -d':' -f2 track_temp.txt > track.txt
    paste output_temp track.txt > results_01.txt

done

sed -i 's/readstats_//g' results_01.txt
sed -i 's/.Read_stats//g' results_01.txt
sed -i 's/\.txt//g' results_01.txt
sed -i 's/\.\///g' results_01.txt

cat results_01.txt

# Transpose the result matrix
cat results_01.txt | Rscript -e 'write.table(t(read.table("stdin",sep="\t",quote="",comment.char="")),sep="\t",quote=F,col.names=F,row.names=F)' > $merged_read_stats_dir/Read_Statistics.${batch_name}.txt

cat $merged_read_stats_dir/Read_Statistics.${batch_name}.txt

rm output_temp
rm track.txt
rm track_temp.txt
rm results_01.txt



cd $current


