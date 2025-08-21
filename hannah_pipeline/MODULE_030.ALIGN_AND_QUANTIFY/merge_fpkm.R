#Run as 
#Rscript /home/uzuny/scripts/merge_fpkm.R sample_dir attribute output_file
args<-commandArgs(TRUE);

sample_dir = args[1]
output_dir = args[2]

#sample_dir = "C:/Users/yuzun/OneDrive - Penn State Health/RESULTS/PROJECTS/2022.IMP.IRINA/FPKM.UNQ/SAMPLES/"
#sample_dir = "C:/Users/yuzun/OneDrive - Penn State Health/RESULTS/PROJECTS/2022.IMP.IRINA/FPKM.UNQ/TEST/"
#output_dir = "C:/Users/yuzun/OneDrive - Penn State Health/RESULTS/PROJECTS/2022.IMP.IRINA/GENE_EXPRESSION_MATRICES/"

#attribute = "gene"

dir.create(output_dir, recursive = T)

print(paste0("Sample directory: ", sample_dir))
print(paste0("Output directory: ", output_dir))

vec_tracking_ids.reference = c()

merge_fpkm_files <- function(samples, attribute)
{
  
  first_sample_name = samples[1]
  
  if(attribute == "gene")
  {
     filename <- paste0(first_sample_name, '/', 'genes.fpkm_tracking')
  }else
  {
    filename <- paste0(first_sample_name, '/', 'isoforms.fpkm_tracking')
  }
  df.fpkm.sample <- read.table(filename, header = T)
  head(df.fpkm.sample)
  vec_tracking_ids.reference = df.fpkm.sample$tracking_id
  df.fpkm.sample$tracking_id[duplicated(df.fpkm.sample$tracking_id)]
  rownames(df.fpkm.sample) = paste0(df.fpkm.sample$tracking_id)
  
  #df.fpkm.merged <- df.fpkm.sample[vec_tracking_ids.reference, 
  #                                 c('tracking_id', 'gene_id', 'gene_short_name', 'locus')]
  #head(df.fpkm.merged)
  
  df.fpkm.annot <- df.fpkm.sample[vec_tracking_ids.reference, 
                                   c('tracking_id', 'gene_id', 'gene_short_name', 'locus')]
  head(df.fpkm.annot)
  
  df.fpkm.annot$tracking_id = substr(df.fpkm.annot$tracking_id, start = 1, stop = 15)
  df.fpkm.annot$gene_id = substr(df.fpkm.annot$gene_id, start = 1, stop = 15)
  
  
  #colnames(merged_df)[5] <- paste0('S_', samples[1])
  list.fpkm = list()
  
  for (i in 1:length(samples) ) {
    sample_name = samples[i]
    print(paste(i, ':' , sample_name)) 
    
    if(attribute == "gene")
    {
      filename <- paste0(sample_name, '/', 'genes.fpkm_tracking')
    }else
    {
      filename <- paste0(sample_name, '/', 'isoforms.fpkm_tracking')
    }
    df.fpkm.sample <- read.table(filename, header = T) 
    head(df.fpkm.sample)

    df.fpkm.sample$tracking_id[duplicated(df.fpkm.sample$tracking_id)]
    rownames(df.fpkm.sample) = paste0(df.fpkm.sample$tracking_id)
    
    df.fpkm.sample = df.fpkm.sample[vec_tracking_ids.reference, c('tracking_id', 'FPKM')]
    print(head(df.fpkm.sample))
    
    vec_fpkm <- as.numeric(df.fpkm.sample$FPKM)
    head(vec_fpkm)
    list.fpkm[[sample_name]] = vec_fpkm
    
    #df.fpkm.sample[df.fpkm.sample$tracking_id == "ENSG00000177301.16", ]
    
    #if(!identical(vec_tracking_ids.reference, vec_tracking_ids.sample))
    #{
    #  stop('Tracking IDs are not matching across samples. Aborting. Check sample FPKM files!!!!!!!!')
    #}
    
    

    
    #df.fpkm.merged = merge(df.fpkm.merged, df.fpkm.sample, by.x = "tracking_id", by.y = "tracking_id")
    
    
    #merged_df <- merge(merged_df, df2, by.x = "tracking_id", by.y = "tracking_id")
  }
  #lapply(list.fpkm, length)
  
  df.fpkms = do.call('cbind', list.fpkm)
  
  head(df.fpkms)
  
  df.fpkm.merged = cbind(df.fpkm.annot, df.fpkms)
  
  head(df.fpkm.merged)
  
  return(df.fpkm.merged)
}


#attribute = "gene"
#setwd('/home/uzuny/gm2/rna_seq/data/fpkm_filtered/all_reads/')
#output_file = paste0("fpkm_combined_expressed.", attribute, ".txt")

setwd(sample_dir)
directories = sort(list.dirs(path = "./") )
print(directories)

samples = c()
for(i in 1:length(directories))
{
  bn = basename(directories[i]) 
  if(bn != ".")
  {
    samples = c(samples, bn)
  }
}
num_samples = length(samples)

print(paste("There are",num_samples,"samples"))

vec_attributes = c('gene', 'transcript')

for(attribute in vec_attributes)
{

  print(attribute)
  print(samples)
  merged_fpkm = merge_fpkm_files(samples, attribute)
  print(head(merged_fpkm))
  print(dim(merged_fpkm))
  
  output_txt = paste0(output_dir, "/FPKM_Matrix.All_",attribute,"s.txt")
  write.table(merged_fpkm, file=output_txt, quote = F, sep='\t', row.names = F, col.names = T)
  output_rds = paste0(output_dir, "/FPKM_Matrix.All_",attribute,"s.rds")
  saveRDS(merged_fpkm, output_rds)
  
  sample_count = length(samples)
  ncol(merged_fpkm)
  merged_fpkm.values = merged_fpkm[,5:(4+sample_count)]
  head(merged_fpkm.values)
  
  temp1 = merged_fpkm.values - rep(1, sample_count)
  temp2 = temp1 > 0
  temp3 = rowSums(temp2)
  temp4 = temp3 > 0
  head(temp4)
  print(paste0("Expressed attributes: ", sum(temp4)))
  merged_fpkm_expressed = merged_fpkm[temp4,] #expressed at least in one sample
  
  output_txt = paste0(output_dir, "/FPKM_Matrix.Expressed_",attribute,"s.txt")
  write.table(merged_fpkm_expressed, file=output_txt, quote = F, sep='\t', row.names = F, col.names = T)
  output_rds = paste0(output_dir, "/FPKM_Matrix.Expressed_",attribute,"s.rds")
  saveRDS(merged_fpkm_expressed, output_rds)

  print(head(merged_fpkm_expressed))
  print(dim(merged_fpkm))
  
}


#process_cell('TH', 'genes')
#process_cell('TR', 'genes')

#setwd('/mnt/isilon/cbmi/tan_lab/uzuny/t1d/rnaseq/fpkm/transcripts.fpkm/')

#process_cell('TH', 'transcripts')
#process_cell('TR', 'transcripts')
