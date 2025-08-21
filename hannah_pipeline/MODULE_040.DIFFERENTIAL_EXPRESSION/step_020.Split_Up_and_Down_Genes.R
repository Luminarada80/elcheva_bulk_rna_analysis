library(openxlsx)

project_name = '2025.NYNRIN.DOD.ELCHEVA'
batch_name = '2025_01'


lab_dir = '/gpfs/Labs/Uzun'

list.files(lab_dir)
  
data_dir = paste0(lab_dir,'/DATA')
metadata_dir = paste0(lab_dir,'/METADATA')
result_dir = paste0(lab_dir,'/RESULTS')
plot_dir = result_dir


gene_annot_dir = paste0(data_dir,'/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT')
biomart_rds =  paste0(gene_annot_dir,"/Biomart_Gene_Table.Clean.rds")

project_data_dir = paste0(data_dir,'/PROJECTS/',project_name)
project_metadata_dir = paste0(metadata_dir,'/PROJECTS/',project_name)

project_result_dir = paste0(result_dir,'/PROJECTS/',project_name)
dea_result_dir = paste0(project_result_dir,'/DIFFERENTIAL_EXPRESSION/',batch_name,'/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest')

gene_type = "coding"
gene_type = "total"
FDR_threshold = 0.05
FC_threshold = 1.5

result_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_',gene_type,'_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.rds')
list_de_results.significant = readRDS(result_rds)

names(list_de_results.significant)

list_de_results.significant$Summary

df_summary = list_de_results.significant$Summary

result_list_new = list()
result_list_new[['Summary']] = df_summary

list_length = length(list_de_results.significant)
count_up_genes = c()
count_down_genes = c()

for(i in 2:list_length)
{
  comparison_name = names(list_de_results.significant)[i]
  df_de = list_de_results.significant[[i]]
  head(df_de)
  print(dim(df_de))
  
  df_de.up = df_de[df_de$log2FC > 0, ]
  head(df_de.up)
  print(nrow(df_de.up))
  count_up_genes[comparison_name] = nrow(df_de.up)
  sheet_name = paste0(comparison_name, '.Up')
  result_list_new[[sheet_name]] = df_de.up
  
  df_de.down = df_de[df_de$log2FC < 0, ]
  head(df_de.down)
  count_down_genes[comparison_name] = nrow(df_de.down)
  sheet_name = paste0(comparison_name, '.Down')
  result_list_new[[sheet_name]] = df_de.down
  
}

df_summary$Sig.DE.Up = count_up_genes
df_summary$Sig.DE.Down = count_down_genes


result_list_new[['Summary']] = df_summary


split_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_', gene_type,
                   '_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.Split.rds')
saveRDS(result_list_new, split_rds)
print(split_rds)

split_xlsx = paste0(dea_result_dir, '/DE_Results.Significantly_DE_', gene_type,
                   '_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.Split.xlsx')
write.xlsx(result_list_new, split_xlsx)
print(split_xlsx)





