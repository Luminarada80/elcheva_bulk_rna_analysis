library("Rsubread")
args <- commandArgs(trailingOnly = T)

#The attribute type which the feature will be built on
quant_level=args[1]
#the path of the annotation file
annotation_file <- args[2]
#the path where you want to save your output
#the path of bam file
bam_file <- args[3]
output_dir <- args[4]
output_name <- args[5]

#bam_file = "/gpfs/Scratch/UzunLab/DATA/PROJECTS/2022.IMP.IRINA/BAM_FILES//shADAR1_INFbeta_1_Aligned.sortedByCoord.out.bam"
#quant_level = "gene"
#annotation_file = "/gpfs/Labs/Uzun/DATA/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/gencode.v38.annotation.gtf"
attrType = "gene_id"

if(quant_level == "gene")
{
  attrType = "gene_id"
}
if(quant_level == "transcript")
{
  attrType = "transcript_id"
}

if(attrType != "gene_id" & attrType != "transcript_id")
{
  stop("Quantification level should be either gene or transcript")
}

temp <-  tail(strsplit(bam_file,"/")[[1]], n=2)
name <- temp[1]
print(bam_file)
print(name)

counts <- featureCounts(bam_file, 
            annot.ext = annotation_file, 
            isGTFAnnotationFile = TRUE,#annotation is in GTF format
            useMetaFeatures = TRUE, #Count reads at the metafeature level
	           GTF.attrType = attrType, #attribute type in the GTF annotation which will 
				     #be used to group features  into meta-features
            isPairedEnd = TRUE,
            nthreads = 1#number of threads for running this function
            )

#output_file_rda = paste(output_dir, output_name ,".rda", sep = "")
#save(counts, file=output_file_rda)

output_file_rds = paste(output_dir, output_name ,".rds", sep = "")
saveRDS(counts, file=output_file_rds)


output_file_txt = paste(output_dir, output_name ,".txt", sep = "")
df = data.frame(rownames(counts$counts), counts$counts)
colnames(df) <- c('id', 'raw_read_count')
head(df)
write.table(df, file=output_file_txt, sep ="\t", quote = F, row.names = F, col.names = T)





