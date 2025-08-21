
library(openxlsx)

project_name='2025.NYNRIN.DOD.ELCHEVA'
batch_name='2025_01'

exp_matrix_dir = paste0("/gpfs/Labs/Uzun/DATA/PROJECTS/",project_name,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name)
project_result_dir = paste0(result_dir,'/PROJECTS/',project_name,'/')

data_type = 'FPKM'
#data_type = 'Read_Count'

exp_matrix_no_annot_gene_symbols_as_rownames_rds = paste0(exp_matrix_dir, '/',data_type,'_Matrix.no_annot.gene_symbols_as_rownames.rds')

df.expression = readRDS(exp_matrix_no_annot_gene_symbols_as_rownames_rds)

gene_means = rowMeans(df.expression)
mean(gene_means)
hist(log10(as.numeric(as.matrix(df.expression))))

df.expression['OAS2', ]
df.expression['OASL', ]
df.expression['IFI44', ]
df.expression['DDX58', ]
df.expression['IRG7', ]



head(df.expression)

cnames = colnames(df.expression)

cnames.interest = cnames
cnames.interest = cnames.interest[grepl('_NT', cnames.interest)]
cnames.interest = cnames.interest[!grepl('ADAR', cnames.interest)]
cnames.interest = cnames.interest[!grepl('Both', cnames.interest)]
cnames.interest


#can_genes = c('OAS2', 'OASL', 'IFI44', 'DDX58', 'IRF7')
can_genes = c('OAS2',  'IFI44', 'DDX58', 'IRF7', 'IRF9', 'STAT1')
df.expression.can = df.expression[can_genes, cnames.interest]
df.expression.can

exp.Ctrl = as.matrix(df.expression.can[, 1:4])
exp.IMP1 = as.matrix(df.expression.can[, 5:8])
exp.IMP2 = as.matrix(df.expression.can[, 9:10])
exp.IMP3 = as.matrix(df.expression.can[, 11:12])

Cntrl_Mean = round(rowMeans(exp.Ctrl), 2)
IMP1_Mean = round(rowMeans(exp.IMP1), 2)
IMP2_Mean = round(rowMeans(exp.IMP2), 2)
IMP3_Mean = round(rowMeans(exp.IMP3), 2)

df.expression.can$Cntrl_Mean = Cntrl_Mean
df.expression.can$IMP1_Mean = IMP1_Mean
df.expression.can$IMP2_Mean = IMP2_Mean
df.expression.can$IMP3_Mean = IMP3_Mean


p.Ctrl_vs_IMP1 =  sapply(1:nrow(exp.Ctrl), function(x) wilcox.test(exp.Ctrl[x,], exp.IMP1[x,])$p.value)
p.Ctrl_vs_IMP2 =  sapply(1:nrow(exp.Ctrl), function(x) wilcox.test(exp.Ctrl[x,], exp.IMP2[x,])$p.value)
p.Ctrl_vs_IMP3 =  sapply(1:nrow(exp.Ctrl), function(x) wilcox.test(exp.Ctrl[x,], exp.IMP3[x,])$p.value)


df.expression.can$p.Cntrl_vs_IMP1 = p.Ctrl_vs_IMP1
df.expression.can$p.Cntrl_vs_IMP2 = p.Ctrl_vs_IMP2
df.expression.can$p.Cntrl_vs_IMP3 = p.Ctrl_vs_IMP3


df.expression.can1 = data.frame(gene = rownames(df.expression.can), df.expression.can)

can.xlsx = paste0(project_result_dir, '/Expression_for_selected_genes.NT.xlsx')
write.xlsx(df.expression.can1, can.xlsx)


par(mfrow = c(2,1))
hist(log10(as.numeric(as.matrix(df.expression[, cnames.interest]))))
hist(log10(as.numeric(as.matrix(df.expression[can_genes, cnames.interest]))))

