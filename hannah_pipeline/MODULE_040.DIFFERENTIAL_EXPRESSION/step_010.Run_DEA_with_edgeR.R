library(openxlsx)
library(edgeR)

project_name='2025.NYNRIN.DOD.ELCHEVA'
batch_name='2025_01'

lab_dir = '/gpfs/Labs/Uzun'

list.files(lab_dir)
  
data_dir = paste0(lab_dir,'/DATA')
metadata_dir = paste0(lab_dir,'/METADATA')
result_dir = paste0(lab_dir,'/RESULTS')
plot_dir = result_dir

gene_annot_dir = paste0(data_dir,'/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/')
biomart_rds =  paste0(gene_annot_dir,"/Biomart_Gene_Table.Clean.rds")

project_data_dir = paste0(data_dir,'/PROJECTS/',project_name)
project_metadata_dir = paste0(metadata_dir,'/PROJECTS/',project_name)

project_result_dir = paste0(result_dir,'/PROJECTS/',project_name,'/')
dir.create(project_result_dir)

exp_matrix_dir = paste0(data_dir,"/PROJECTS/",project_name,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name,"/")

comparison_list_txt = paste0(project_metadata_dir,'/BULK_RNA_SEQ/BATCH_2024_04/Comparisons.BATCH_',batch_name,'.txt') 

sample_annotation_txt = paste0(project_metadata_dir,'/BULK_RNA_SEQ/BATCH_2024_04/Sample_Annotation.BATCH_',batch_name'.txt') 
fpkm.rds = paste0(project_data_dir,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name,"/FPKM_Matrix.All_genes.rds")
read_count.rds = paste0(project_data_dir,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name,"/Read_count_matrix.All_genes.rds")

df.comparisons = read.table(comparison_list_txt, sep = '\t', header = T)

df.sample_annotation = read.table(sample_annotation_txt, sep = '\t', header = T)
rownames(df.sample_annotation) = df.sample_annotation$Sample_Name
head(df.sample_annotation)

df.biomart = readRDS(biomart_rds)
head(df.biomart)
dim(df.biomart)

df.biomart[grepl('ENSG00000102003', df.biomart$ensembl_gene_id_version), ]

df.fpkm = readRDS(fpkm.rds)
df.fpkm[1:5, 1:5]
dim(df.fpkm)
df.fpkm[duplicated(df.fpkm$tracking_id), 'tracking_id']

numcol = ncol(df.fpkm)
df.fpkm[,5:16] = round(df.fpkm[,5:16], 2)

temp = df.fpkm$gene_short_name
temp[grepl('RIG', temp)]


colnames(df.fpkm)[1] = "Gene_id"
df.fpkm[1:5, 1:5]
#df.fpkm$Gene_id = NULL
#df.fpkm[1:5, 1:5]

df.read_counts = readRDS(read_count.rds)
#read_counts_txt = paste0(project_data_dir, '/Gene_Expression_Matrix_Read_Counts.txt')
#df.read_counts = read.table(read_counts_txt, header = T, sep = '\t')

df.read_counts[1:5, 1:5]
dim(df.read_counts)

colnames(df.read_counts)


dim(df.read_counts)


n_comp = nrow(df.comparisons)
n_annot = nrow(df.sample_annotation)
head(df.sample_annotation)

Sample_Names = df.sample_annotation$Sample_Name

test_type = 'glmQLFTest'

#gene_type = 'all_gene_types'
#gene_type = 'protein_coding'

#filtering = 'None'
filtering = 'FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP'

dea_result_dir = paste0(project_result_dir, '/DIFFERENTIAL_EXPRESSION/', batch_name,'/', filtering , '/',test_type, '/')
dir.create(dea_result_dir, recursive = T)



list_de_results = list()
list_de_results.protein_coding = list()

#count_DE.all = c(rep(0, n_comp))
#count_DE.protein_coding = c(rep(0, n_comp))

count_DE.all = c()
count_DE.protein_coding = c()
count_expressed.all = c()
count_expressed.protein_coding = c()




list_de_results[['Summary']] = df.comparisons
list_de_results.protein_coding[['Summary']] = df.comparisons

#n_comp = 1
#edger uses log2 FC. 
FC_threshold = 1.5
logFC_threshold = log2(FC_threshold) #%50 higher expression

FDR_threshold = 0.05

for(i in 1:n_comp)
{
  #print(i)
  comparison_id =  df.comparisons[i, 'Comparison_Id']
  print(comparison_id)
  
  group_name_1 = as.character(df.comparisons[i, 'Group_1'])
  group_name_2 = as.character(df.comparisons[i, 'Group_2'])
  

  comparison_name = paste0(group_name_1, '_vs_', group_name_2)
  print(paste(comparison_id, ':' , comparison_name))
  
  sample_groups = as.character(df.sample_annotation$Group_Name)

  group_1_bool = (sample_groups == group_name_1 )
  group_1_sample_names = Sample_Names[group_1_bool]  
  
  group_2_bool = (sample_groups == group_name_2 )
  group_2_sample_names = Sample_Names[group_2_bool]  
  
  #filter genes
  

  df.fpkm.comparison = df.fpkm[, c(group_1_sample_names, group_2_sample_names)]
  head(df.fpkm.comparison)
  df.fpkm.comparison.two_digits = round(df.fpkm.comparison, 2)
  head(df.fpkm.comparison.two_digits)
  
  
  df.fpkm.g1 = df.fpkm[, group_1_sample_names]
  df.fpkm.g2 = df.fpkm[, group_2_sample_names]
  
  df.fpkm.g1 = round(df.fpkm.g1, 2)
  df.fpkm.g2 = round(df.fpkm.g2, 2)
  
  head(df.fpkm.g1)
  head(df.fpkm.g2)
  
  g1_fpkm_mean = round(rowMeans(df.fpkm.g1), 2)
  g2_fpkm_mean = round(rowMeans(df.fpkm.g2) ,2)
  head(g1_fpkm_mean)
  head(g2_fpkm_mean)
  
  
  bool.fpkm_gt_one.g1 = df.fpkm.g1 > 1 
  bool.fpkm_gt_one.g2 = df.fpkm.g2 > 1 
  
  head(bool.fpkm_gt_one.g1)
  head(bool.fpkm_gt_one.g2)
  
  count.fpkm_gt_one.g1 = rowSums(bool.fpkm_gt_one.g1)
  count.fpkm_gt_one.g2 = rowSums(bool.fpkm_gt_one.g2)
  
  head(count.fpkm_gt_one.g1)
  head(count.fpkm_gt_one.g2)
  
  bool.fpkm_gt_one_at_least_two_samples.g1 = count.fpkm_gt_one.g1 > 1
  bool.fpkm_gt_one_at_least_two_samples.g2 = count.fpkm_gt_one.g2 > 1
  
  head(bool.fpkm_gt_one_at_least_two_samples.g1)
  head(bool.fpkm_gt_one_at_least_two_samples.g2)
  
  bool.fpkm_gt_one_at_least_two_samples_in_either_groups = bool.fpkm_gt_one_at_least_two_samples.g1 | bool.fpkm_gt_one_at_least_two_samples.g2
  gene_id_and_version.expressed = rownames(df.fpkm.comparison[bool.fpkm_gt_one_at_least_two_samples_in_either_groups, ])
  
  head(gene_id_and_version.expressed)
  length(gene_id_and_version.expressed)
  
  
  df.read_counts.expressed = df.read_counts[gene_id_and_version.expressed, ]
  dim(df.read_counts.expressed)
  #select samples
  
  if(filtering == 'FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP')
  {
    df.read_counts.comparison = df.read_counts.expressed[, c(group_1_sample_names, group_2_sample_names)]
  }else
  {
    df.read_counts.comparison = df.read_counts.expressed[, c(group_1_sample_names, group_2_sample_names)]
  }
  
  head(df.read_counts.comparison)
  dim(df.read_counts.comparison)
  
  sum(is.na(df.read_counts.comparison))
  


  
  #build the model
  group = c(rep("G1", length(group_1_sample_names)), rep("G2", length(group_2_sample_names)))
  batch = df.sample_annotation[c(group_1_sample_names, group_2_sample_names), 'Batch_ID']
 
 
  head(df.read_counts.comparison)
  
  
  #apply TMM normalization
  y <- DGEList(counts = df.read_counts.comparison, group = group )
  y <- calcNormFactors(y)	
  
  if(test_type == 'exactTest')
  {
    #Estimate the dispersion, necessary#
    y <- estimateCommonDisp(y, verbose=T)
    y <- estimateTagwiseDisp(y)
    
    #Differential expression
    et <- exactTest(y)
    dim(et)
    head(et)
    test_result = et
  }
  if(test_type == 'glmQLFTest')
  {
    design_no_batch <- model.matrix(~group)
    
    if(length(batch) > 0)
    {
      # There are multiple batches
      design_with_batch <- model.matrix(~group+batch)
      
      y <- estimateDisp(y, design_with_batch)
      
      #merged_vars <- factor(paste(group, batch, sep="."))
      
      fit <- tryCatch({
            fitx = glmQLFit(y, design_with_batch)
            fitx
          }, warning = function(w) {
            fitx = glmQLFit(y,design_with_batch)
            print(w)
            fitx
          }, error = function(e) {
            y <- estimateDisp(y,design_no_batch)
            print(e)
            fitx = glmQLFit(y,design_no_batch)
            
            fitx
          }, finally = {
            
          }
        )
      
    }else{
      #All samples are from the same batch or there is no batch information
      y <- estimateDisp(y, design_no_batch)
      fit = glmQLFit(y, design_no_batch)
    }
    #fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit, coef=2)
    test_result = qlf
    
  }
  
  n_genes = nrow(test_result)
  top <- topTags(test_result, n=n_genes)
  
  
  head(top)
  dim(top)
  
  #df.fpkm.comparison.two_digits['ENSG00000089199.10', ]
  #de <- subset(top$table,top$table$FDR <= 1)
  #head(de)
  top$adjust.method
  top$test
  df_de_result = top$table
  head(df_de_result)
  dim(df_de_result)
  
  df_de_result$logFC = round(df_de_result$logFC, 2)
  df_de_result$logCPM = round(df_de_result$logCPM, 2)
  #df_de_result$PValue = format(df_de_result$PValue, scientific = T, digits = 3)
  #df_de_result$FDR = format(df_de_result$FDR, scientific = T, digits = 3)
  
  df_de_result.1 = data.frame(gene_id = rownames(df_de_result), df_de_result)
  head(df_de_result.1)
  dim(df_de_result.1)
  sum(df_de_result.1$FDR < 0.01)
  sum(df_de_result.1$FDR < 0.05)
  sum(df_de_result.1$FDR < 0.10)
  
  sum(df_de_result.1$FDR < 0.01 & abs(df_de_result.1$logFC) > 0.5)
  sum(df_de_result.1$FDR < 0.01 & abs(df_de_result.1$logFC) > 1)
  
  
  count_de_0.05 = sum(df_de_result$FDR < 0.05)

  msg = paste('FDR 0.05:', count_de_0.05) 
  print(msg)
  
  sum(df_de_result$PValue < 0.01)
  sum(df_de_result$PValue < 0.05)
  sum(df_de_result$PValue < 0.10)
  
  gene_ids = as.character(df_de_result.1$gene_id)
  
  gene_ids.wo_version = gsub("\\..*","",gene_ids)
  head(gene_ids.wo_version)
  head(df.biomart)
  df.gene_annot = df.biomart[gene_ids.wo_version, c('gene_symbol', 'gene_biotype', 'gene_desc', 'chrom_name')]
  head(df.gene_annot)
  df.gene_annot[df.gene_annot$gene_symbol == 'RIGI', ]
  #df.gene_annot['ENSG00000102003',  ]
  
  
  df_de_result.2 = data.frame(df.gene_annot, df_de_result.1,
                              G1_MEAN_FPKM = g1_fpkm_mean[gene_ids], 
                              G2_MEAN_FPKM = g2_fpkm_mean[gene_ids], 
                              df.fpkm.comparison[gene_ids, ])
  head(df_de_result.2)
  
  all_p_values = df_de_result.1$PValue
  names(all_p_values) = df_de_result.1$gene_id
  head(all_p_values)
  
  gene_ids.all = as.character(df_de_result.2$gene_id)
  

  
  df_de_result.p_value_sorted = df_de_result.2[order(df_de_result.2$PValue), ]
  
  dim(df_de_result.p_value_sorted)
  length(unique(df_de_result.p_value_sorted$gene_symbol))
  length(unique(df_de_result.p_value_sorted$gene_id))
  
  de_result.unique = df_de_result.p_value_sorted[!duplicated(df_de_result.p_value_sorted$gene_symbol), ] 
  count_expressed.all[comparison_name] = nrow(de_result.unique)
  list_de_results[[comparison_name]] = de_result.unique

  de_count.all = sum(de_result.unique$FDR < FDR_threshold & abs(de_result.unique$logFC) > logFC_threshold)
  count_DE.all[comparison_name] = de_count.all
  
  
  
  
  
  #Filter coding genes
  df_de_result.4 = df_de_result.2
  df_de_result.4$gene_biotype = as.character(df_de_result.4$gene_biotype)
  df_de_result.4$gene_biotype[is.na(df_de_result.4$gene_biotype)] = 'not_available'
  df_de_result.protein_coding =df_de_result.4[df_de_result.4$gene_biotype == 'protein_coding', ]
  
  head(df_de_result.protein_coding)
  

  df_de_result.protein_coding.p_value_sorted = df_de_result.protein_coding[order(df_de_result.protein_coding$PValue), ]
  de_result_coding.unique = df_de_result.protein_coding.p_value_sorted[!duplicated(df_de_result.protein_coding.p_value_sorted$gene_symbol), ] 
  count_expressed.protein_coding[comparison_name] = nrow(de_result_coding.unique)
  list_de_results.protein_coding[[comparison_name]] = de_result_coding.unique

  de_count.protein_coding = sum(de_result_coding.unique$FDR < FDR_threshold & abs(de_result_coding.unique$logFC) > logFC_threshold)
  count_DE.protein_coding[comparison_name] = de_count.protein_coding
  
}


df_de_counts.all = cbind(df.comparisons, 'Expressed.All' = count_expressed.all,  'Sig_DE.All' = count_DE.all)
df_de_counts.protein_coding = cbind(df.comparisons, 'Expressed.Protein_Coding' = count_expressed.protein_coding, 'Sig_DE.Protein_Coding' = count_DE.protein_coding)

list_de_results[['Summary']] = df_de_counts.all
list_de_results.protein_coding[['Summary']] = df_de_counts.protein_coding



result_rds = paste0(dea_result_dir, '/DE_Results.All_expressed_coding_genes.rds')
saveRDS(list_de_results.protein_coding, result_rds)
#list_de_results.protein_coding = readRDS(result_rds)

result_rds = paste0(dea_result_dir, '/DE_Results.All_expressed_total_genes.rds')
saveRDS(list_de_results, result_rds)
#list_de_results = readRDS(result_rds)



list_de_results.significant = list()
list_de_results.significant.protein_coding = list()

list_de_results.significant[['Summary']] = list_de_results[['Summary']]
list_de_results.significant.protein_coding[['Summary']] = list_de_results.protein_coding[['Summary']]


for(i in 2:length(list_de_results))
{
  #All expressed genes
  df.de = list_de_results[[i]]
  head(df.de)
  df.de$logCPM = NULL
  df.de$F = NULL
  cnames = colnames(df.de)
  cnames = gsub('logFC', 'log2FC', cnames)
  colnames(df.de) = cnames
  list_de_results[[i]] = df.de
  df.de.sig = df.de[df.de$FDR < FDR_threshold & abs(df.de$log2FC) > logFC_threshold, ]
  dim(df.de.sig) 
  list_de_results.significant[[i]] = df.de.sig
  
  
  #Coding genes
  df.de = list_de_results.protein_coding[[i]]
  head(df.de)
  df.de$logCPM = NULL
  df.de$F = NULL
  cnames = colnames(df.de)
  cnames = gsub('logFC', 'log2FC', cnames)
  colnames(df.de) = cnames
  list_de_results.protein_coding[[i]] = df.de
  df.de.sig = df.de[df.de$FDR < FDR_threshold & abs(df.de$log2FC) > logFC_threshold, ]
  dim(df.de.sig) 
  list_de_results.significant.protein_coding[[i]] = df.de.sig
  
}


sheet_names = names(list_de_results)
sheet_names= gsub('Cntrl', 'CNT', sheet_names)
sheet_names= gsub('__vs__', '.vs.', sheet_names)
sheet_names= gsub('beta', 'B', sheet_names)

names(list_de_results) = sheet_names
names(list_de_results.protein_coding) = sheet_names
names(list_de_results.significant) = sheet_names
names(list_de_results.significant.protein_coding) = sheet_names


result_xlsx = paste0(dea_result_dir, '/DE_Results.All_expressed_coding_genes.xlsx')
write.xlsx(list_de_results.protein_coding, result_xlsx)

result_xlsx = paste0(dea_result_dir, '/DE_Results.All_expressed_total_genes.xlsx')
write.xlsx(list_de_results, result_xlsx)


result_xlsx = paste0(dea_result_dir, '/DE_Results.Significantly_DE_coding_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.xlsx')
write.xlsx(list_de_results.significant.protein_coding, result_xlsx)

result_xlsx = paste0(dea_result_dir, '/DE_Results.Significantly_DE_total_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.xlsx')
write.xlsx(list_de_results.significant, result_xlsx)



result_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_coding_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.rds')
saveRDS(list_de_results.significant.protein_coding, result_rds)

result_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_total_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.rds')
saveRDS(list_de_results.significant, result_rds)













