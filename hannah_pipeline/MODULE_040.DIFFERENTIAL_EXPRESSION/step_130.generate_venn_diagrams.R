library(venn)
library(openxlsx)

round_no = '1_and_2'


research_dir = '/gpfs/Labs/Uzun/'

list.files(research_dir)
  
data_dir = paste0(research_dir, '/DATA/')
metadata_dir = paste0(research_dir, '/METADATA/')
result_dir = paste0(research_dir, '/RESULTS/')
plot_dir = result_dir


project_name = '/2022.IMP.IRINA/'


project_data_dir = paste0(data_dir, '/PROJECTS/', project_name, '/')
project_metadata_dir = paste0(metadata_dir, '/PROJECTS/',project_name, '/')

project_result_dir = paste0(result_dir, '/PROJECTS/',  project_name, '/')
dir.create(project_result_dir)

plot_dir = paste0(project_result_dir, '/PLOTS/')
dir.create(plot_dir)

comparison_list_txt = paste0(project_metadata_dir, '/Comparisons.Round_',round_no,'.txt') 
df.comparisons = read.table(comparison_list_txt, sep = '\t', header = T)
rownames(df.comparisons) = df.comparisons$Comparison_Id

df.comparisons$comparison_name = paste0(df.comparisons$G1_Genotype, '_' ,
                                        df.comparisons$G1_Condition, '.vs.', 
                                        df.comparisons$G2_Genotype, '_' , 
                                        df.comparisons$G2_Condition)

df.comparisons$comparison_name= gsub('Cntrl', 'CNT', df.comparisons$comparison_name)
df.comparisons$comparison_name= gsub('__vs__', '.vs.', df.comparisons$comparison_name)
df.comparisons$comparison_name= gsub('beta', 'B', df.comparisons$comparison_name)

cross_list_txt = paste0(project_metadata_dir, '/Cross_Comparisons.txt') 
df.cross = read.table(cross_list_txt, sep = '\t', header = T)

gene_biotypes = 'total'
gene_biotypes = 'coding'


test_type = 'glmQLFTest'

#gene_type = 'all_gene_types'
#gene_type = 'protein_coding'

#filtering = 'None'
filtering = 'FPKM_gt_one_in_at_least_two_samples_in_either_group'


plot_dir = paste0(dea_result_dir, '/VENN.', toupper(gene_biotypes), '/')
dir.create(plot_dir)

FC_threshold = 1.5
log2FC_threshold = log2(FC_threshold) 

FDR_threshold = 0.05

dea_result_dir = paste0(project_result_dir, '/DEA.ROUND_', round_no,'/', filtering , '/',test_type, '/')
result_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_',gene_biotypes,'_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.rds')
list.sig_de = readRDS(result_rds)
n_cross = nrow(df.cross)

list.cross = list()

#direction = 'All_sig'
#direction = 'Up'
#direction = 'Down'


directions = c('All_sig', 'Up', 'Down')

for(direction in directions)
{
  pdf_file = paste0(plot_dir,'DEG.Venn.', direction,'.pdf')
  print(pdf_file)
  pdf(pdf_file, width = 12, height = 4)
  par(mfrow = c(1,3))
  
  for(i in 1:n_cross)
  {
    print(i)
    cross_id = df.cross$Cross_Id[i]
    comp.A = df.cross$Set_A[i]
    comp.B = df.cross$Set_B[i]
    
    cross_name = paste0(comp.A, '_x_', comp.B)
    print(cross_name)
  
    comparison_name_A = df.comparisons[comp.A, 'comparison_name'] 
    comparison_name_B = df.comparisons[comp.B, 'comparison_name'] 
    
    df.de.A = list.sig_de[[comparison_name_A]]
    df.de.B = list.sig_de[[comparison_name_B]]
    
    df.de.A[df.de.A$gene_symbol == 'OAS2', ]
    
    if(direction == 'Up')
    {
      df.de.A = df.de.A[df.de.A$log2FC > 0, ]
      df.de.B = df.de.B[df.de.B$log2FC > 0, ]
    }
    
    if(direction == 'Down')
    {
      df.de.A = df.de.A[df.de.A$log2FC < 0, ]
      df.de.B = df.de.B[df.de.B$log2FC < 0, ]    
    }
    
    head(df.de.A)
    head(df.de.B)
    
    gene_ids.A = rownames(df.de.A)
    gene_ids.B = rownames(df.de.B)
    
    head(gene_ids.A)
    head(gene_ids.B)
    
    gene_ids.shared = intersect(gene_ids.A, gene_ids.B)
    
    head(gene_ids.shared)
    length(gene_ids.shared)
    
    df.de.A.shared = df.de.A[gene_ids.shared, ]
    cnames <- colnames(df.de.A.shared)
    cnames[6:10] = paste0(comp.A, '.', cnames[6:10])
    cnames -> colnames(df.de.A.shared)
    
    df.de.B.shared = df.de.B[gene_ids.shared, ]
    cnames <- colnames(df.de.B.shared)
    cnames[6:10] = paste0(comp.B, '.', cnames[6:10])
    cnames -> colnames(df.de.B.shared)
    
    #df.de.B.shared[, 1:5] = NULL
    
    df.de.shared = cbind(df.de.A.shared[, 1:10], df.de.B.shared[, 6:10])
    
    df.de.shared$min_p = pmin(df.de.shared[, 7], df.de.shared[, 12])
    head(df.de.shared)
    df.de.shared = df.de.shared[order(df.de.shared$min_p), ]
    df.de.shared$min_p = NULL
    
    list.cross[[cross_name]] = df.de.shared
    set_names = c(comparison_name_A, comparison_name_B)
    set_names = gsub('shCNT_NT.vs.', '', set_names)
    
  
    
    venn(x = list(gene_ids.A, gene_ids.B), 
         snames = set_names, ilcs = 2, sncs = 2, box = F,
         zcolor = "style")
    
    labels = paste0('shCNT_NT vs :')
    mtext(text = labels, line = -2, cex = 1)
    
  }
  
  dev.off()
  
  
  
  rds_file = paste0(dea_result_dir, '/Shared_DE.',gene_biotypes,'_genes.FDR_',FDR_threshold, '.FC_',FC_threshold, '.', direction ,'.rds')
  saveRDS(list.cross, rds_file)
  
  xlsx_file = paste0(dea_result_dir, '/Shared_DE.',gene_biotypes,'_genes.FDR_',FDR_threshold, '.FC_',FC_threshold, '.', direction ,'.xlsx')
  write.xlsx(list.cross, xlsx_file)
  

}








