rm(list = ls())

library(data.table)

project_name='2025.NYNRIN.DOD.ELCHEVA'
batch_name='2025_01'

exp_matrix_dir = paste0("/gpfs/Labs/Uzun/DATA/PROJECTS/",project_name,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name)

data_type = 'FPKM'
#data_type = 'Read_Count'

exp_matrix_rds = paste0(exp_matrix_dir,'/',data_type,'_Matrix.All_genes.rds')

exp_matrix_with_annot_rds = paste0(exp_matrix_dir,'/',data_type,'_Matrix.cleaned.with_annot.rds')
exp_matrix_no_annot_rds = paste0(exp_matrix_dir,'/',data_type,'_Matrix.cleaned.no_annot.rds')
exp_matrix_no_annot_gene_symbols_as_rownames_rds = paste0(exp_matrix_dir,'/',data_type,'_Matrix.cleaned.no_annot.gene_symbols_as_rownames.rds')

df.expression = readRDS(exp_matrix_rds)
head(df.expression)
dim(df.expression)
length(unique(df.expression$gene_id))

df.expression[duplicated(df.expression$gene_short_name), ]


df.expression[duplicated(df.expression$gene_id), 1:10]
bool.PAR_Y_ids = grepl('PAR_Y', rownames(df.expression)) 

df.expression[bool.PAR_Y_ids, 1:20]
df.expression[bool.PAR_Y_ids, 21:ncol(df.expression)]


df.expression.PAR_Y_removed = df.expression[!bool.PAR_Y_ids,]

df.expression.PAR_Y_removed[duplicated(df.expression.PAR_Y_removed$gene_id), 1:10]

head(df.expression.PAR_Y_removed)

saveRDS(df.expression.PAR_Y_removed, exp_matrix_with_annot_rds)



all_cols = colnames(df.expression.PAR_Y_removed)
annot_cols = c( 'tracking_id', 'gene_id', 'gene_short_name',
                'locus')
value_cols = sort(setdiff(all_cols, annot_cols))


df.expression.values = df.expression.PAR_Y_removed[, value_cols]
if(data_type == 'FPKM')
{
  df.expression.values = round(df.expression.values, 2)
}

head(df.expression.values)

saveRDS(df.expression.values, exp_matrix_no_annot_rds)


mean_expression = rowMeans(df.expression[, value_cols])
dim(df.expression)
head(mean_expression)
length(mean_expression)
df.expression.sorted_by_mean_exp = df.expression[order(mean_expression, decreasing = T), ]
head(df.expression.sorted_by_mean_exp)
tail(df.expression.sorted_by_mean_exp)

df.expression.unique_gene_symbols = df.expression.sorted_by_mean_exp[!duplicated(df.expression.sorted_by_mean_exp$gene_short_name), ] 
dim(df.expression.unique_gene_symbols)
rownames(df.expression.unique_gene_symbols) = df.expression.unique_gene_symbols$gene_short_name
head(df.expression.unique_gene_symbols)

df.expression.no_annot.gene_symbols_as_rownames = df.expression.unique_gene_symbols[, value_cols]
head(df.expression.no_annot.gene_symbols_as_rownames)

saveRDS(df.expression.unique_gene_symbols, exp_matrix_no_annot_gene_symbols_as_rownames_rds)


