library(openxlsx)
library(edgeR)

project_name='2025.NYNRIN.DOD.ELCHEVA'
batch_name='2025_01'

lab_dir = '/gpfs/Labs/Uzun'

data_dir      = paste0(lab_dir,'/DATA')
metadata_dir  = paste0(lab_dir,'/METADATA')
result_dir    = paste0(lab_dir,'/RESULTS')
plot_dir      = result_dir

gene_annot_dir = paste0(data_dir,'/GENOMES/ANNOTATION/HUMAN/HG38/GENE_ANNOT/')
biomart_rds    = paste0(gene_annot_dir,"/Biomart_Gene_Table.Clean.rds")

project_data_dir     = paste0(data_dir,'/GRANT_APPS/',project_name)
project_metadata_dir = paste0(metadata_dir,'/GRANT_APPS/',project_name)
project_result_dir   = paste0(result_dir,'/GRANT_APPS/',project_name,'/')

exp_matrix_dir = paste0(data_dir,"/GRANT_APPS/",project_name,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name,"/")

comparison_list_txt  = paste0(project_metadata_dir,'/Comparisons.BATCH_',batch_name,'.txt') 
sample_annotation_txt= paste0(project_metadata_dir,'/Sample_Annotation.BATCH_',batch_name,'.txt') 

### CHANGED: we will not use FPKM matrices; we operate on counts + TMM
# fpkm.rds     = paste0(project_data_dir,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name,"/FPKM_Matrix.All_genes.rds")
read_count.rds = paste0(project_data_dir,"/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/",batch_name,"/Read_Count_Matrix.All_genes.rds")

df.comparisons      = read.table(comparison_list_txt, sep = '\t', header = TRUE, check.names = FALSE)
df.sample_annotation= read.table(sample_annotation_txt, sep = '\t', header = TRUE, check.names = FALSE)
rownames(df.sample_annotation) = df.sample_annotation$Sample_Name

df.biomart = readRDS(biomart_rds)  # rownames should be Ensembl IDs (no version) or has ensembl_gene_id_version

### Load counts (genes x samples). Keep raw counts for edgeR.
df.read_counts = readRDS(read_count.rds)

### Harmonize sample names (strip .gene and _S##) to match Sample_Annotation
normalize_id <- function(x) {
  x <- sub("\\.gene$", "", x)
  x <- sub("_S\\d+$", "", x)
  x
}
colnames(df.read_counts) <- normalize_id(colnames(df.read_counts))

### Keep only annotated samples (and in defined order)
Sample_Names <- df.sample_annotation$Sample_Name
missing_c <- setdiff(Sample_Names, colnames(df.read_counts))
if (length(missing_c)) stop("Missing in counts: ", paste(missing_c, collapse = ", "))
df.read_counts <- df.read_counts[, Sample_Names, drop = FALSE]

### Strip Ensembl versions from rownames to join with biomart later
rownames(df.read_counts) <- sub("\\..*$", "", rownames(df.read_counts))

n_comp = nrow(df.comparisons)

test_type = 'glmQLFTest'
FC_threshold  = 1.5
logFC_threshold = log2(FC_threshold)
FDR_threshold = 0.05

dea_result_dir = paste0(project_result_dir, '/DIFFERENTIAL_EXPRESSION/', batch_name,'/',
                        'CPM_FILTER_TMM', '/', test_type, '/')
dir.create(dea_result_dir, recursive = TRUE, showWarnings = FALSE)

list_de_results                       = list()
list_de_results.protein_coding        = list()
list_de_results[['Summary']]          = df.comparisons
list_de_results.protein_coding[['Summary']] = df.comparisons

count_DE.all              = c()
count_DE.protein_coding   = c()
count_expressed.all       = c()
count_expressed.protein_coding = c()

### Helper: pull biomart rows by Ensembl (no version), with safe fallback
bm_pick <- function(ids_novers) {
  ids_novers <- as.character(ids_novers)
  # If df.biomart has a column 'ensembl_gene_id_version', derive no-version index
  if (!is.null(df.biomart$ensembl_gene_id_version)) {
    idx <- sub("\\..*$", "", df.biomart$ensembl_gene_id_version)
    bm2 <- df.biomart
    rownames(bm2) <- idx
    return(bm2[ids_novers, c('gene_symbol','gene_biotype','gene_desc','chrom_name')])
  } else {
    # assume rownames are ensembl IDs (no version)
    return(df.biomart[ids_novers, c('gene_symbol','gene_biotype','gene_desc','chrom_name')])
  }
}

for (i in 1:n_comp) {
  comparison_id =  df.comparisons[i, 'Comparison_Id']
  group_name_1  = as.character(df.comparisons[i, 'Group_1'])
  group_name_2  = as.character(df.comparisons[i, 'Group_2'])
  comparison_name = paste0(group_name_1, '_vs_', group_name_2)
  message(paste(comparison_id, ":", comparison_name))

  sample_groups <- as.character(df.sample_annotation$Group_Name)
  group_1_sample_names <- Sample_Names[sample_groups == group_name_1]
  group_2_sample_names <- Sample_Names[sample_groups == group_name_2]

  ## ---- Build comparison count matrix (genes x samples)
  df.read_counts.comparison <- df.read_counts[, c(group_1_sample_names, group_2_sample_names), drop = FALSE]

  ## ---- edgeR: DGEList + TMM normalization
  group <- c(rep("G1", length(group_1_sample_names)), rep("G2", length(group_2_sample_names)))
  batch <- df.sample_annotation[c(group_1_sample_names, group_2_sample_names), 'Batch_ID']
  group <- factor(group)

  y <- DGEList(counts = df.read_counts.comparison, group = group)
  y <- calcNormFactors(y, method = "TMM")     # library-size normalization

  ## ---- CHANGED: Filter by TMM-CPM, not FPKM
  cpm_tmm <- cpm(y, normalized.lib.sizes = TRUE)  # genes x samples
  keep_g1 <- rowSums(cpm_tmm[, group == "G1", drop = FALSE] >= 1) >= 2
  keep_g2 <- rowSums(cpm_tmm[, group == "G2", drop = FALSE] >= 2)
  keep    <- keep_g1 | keep_g2

  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")     # re-normalize after filtering

  # Store post-filter CPMs and group means for reporting
  cpm_tmm_post <- cpm(y, normalized.lib.sizes = TRUE)  # genes x samples
  g1_cpm_mean  <- rowMeans(cpm_tmm_post[, group == "G1", drop = FALSE])
  g2_cpm_mean  <- rowMeans(cpm_tmm_post[, group == "G2", drop = FALSE])

  ## ---- Dispersion + QL fit (with/without batch)
  design_no_batch   <- model.matrix(~ group)
  use_batch <- !all(is.na(batch)) && length(unique(batch)) > 1
  if (use_batch) {
    batch <- factor(batch)
    design_with_batch <- model.matrix(~ group + batch)
    y  <- estimateDisp(y, design_with_batch)
    fit <- glmQLFit(y, design_with_batch)
    qlf <- glmQLFTest(fit, coef = 2)  # groupG2 vs groupG1
  } else {
    y  <- estimateDisp(y, design_no_batch)
    fit <- glmQLFit(y, design_no_batch)
    qlf <- glmQLFTest(fit, coef = 2)
  }

  top <- topTags(qlf, n = Inf)
  df_de_result <- top$table
  df_de_result$logFC  <- round(df_de_result$logFC,  2)
  df_de_result$logCPM <- round(df_de_result$logCPM, 2)

  df_de_result.1 <- data.frame(gene_id = rownames(df_de_result), df_de_result,
                               row.names = NULL, check.names = FALSE)

  ## ---- Annotation: map Ensembl (no version) -> biomart
  gene_ids            <- df_de_result.1$gene_id
  gene_ids_novers     <- sub("\\..*$", "", gene_ids)
  df.gene_annot       <- bm_pick(gene_ids_novers)

  ## ---- CHANGED: replace mean FPKM columns with mean TMM-CPM
  # Also bind per-sample CPMs if you want (commented)
  per_sample_cpm <- t(cpm_tmm_post[rownames(df_de_result), , drop = FALSE]) # samples x genes (optional view)
  df_de_result.2 <- data.frame(
    df.gene_annot,
    df_de_result.1,
    G1_MEAN_CPM_TMM = g1_cpm_mean[rownames(df_de_result)],
    G2_MEAN_CPM_TMM = g2_cpm_mean[rownames(df_de_result)]
    # , t(cpm_tmm_post[rownames(df_de_result), , drop = FALSE])  # uncomment to append per-sample CPMs
  )
  rownames(df_de_result.2) <- NULL

  ## ---- Unique by gene symbol (like your original)
  df_de_result.2$gene_symbol <- as.character(df_de_result.2$gene_symbol)
  df_de_result.2$gene_biotype <- as.character(df_de_result.2$gene_biotype)
  df_de_result.p_value_sorted <- df_de_result.2[order(df_de_result.2$PValue), ]
  de_result.unique <- df_de_result.p_value_sorted[!duplicated(df_de_result.p_value_sorted$gene_symbol), ]

  count_expressed.all[comparison_name] <- nrow(de_result.unique)
  list_de_results[[comparison_name]]   <- de_result.unique

  de_count.all <- sum(de_result.unique$FDR < FDR_threshold & abs(de_result.unique$logFC) > logFC_threshold)
  count_DE.all[comparison_name] <- de_count.all

  ## ---- Protein-coding subset
  df_de_result.4 <- df_de_result.2
  df_de_result.4$gene_biotype[is.na(df_de_result.4$gene_biotype)] <- 'not_available'
  df_de_result.protein_coding <- df_de_result.4[df_de_result.4$gene_biotype == 'protein_coding', , drop = FALSE]
  df_de_result.protein_coding.p_value_sorted <- df_de_result.protein_coding[order(df_de_result.protein_coding$PValue), ]
  de_result_coding.unique <- df_de_result.protein_coding.p_value_sorted[!duplicated(df_de_result.protein_coding.p_value_sorted$gene_symbol), ]

  count_expressed.protein_coding[comparison_name] <- nrow(de_result_coding.unique)
  list_de_results.protein_coding[[comparison_name]] <- de_result_coding.unique

  de_count.protein_coding <- sum(de_result_coding.unique$FDR < FDR_threshold & abs(de_result_coding.unique$logFC) > logFC_threshold)
  count_DE.protein_coding[comparison_name] <- de_count.protein_coding
}

## ---- Summary tabs
df_de_counts.all = cbind(df.comparisons, 'Expressed.All' = count_expressed.all,  'Sig_DE.All' = count_DE.all)
df_de_counts.protein_coding = cbind(df.comparisons, 'Expressed.Protein_Coding' = count_expressed.protein_coding,
                                    'Sig_DE.Protein_Coding' = count_DE.protein_coding)

list_de_results[['Summary']] <- df_de_counts.all
list_de_results.protein_coding[['Summary']] <- df_de_counts.protein_coding

## ---- Significant-only lists (preserve your structure)
list_de_results.significant <- list()
list_de_results.significant.protein_coding <- list()
list_de_results.significant[['Summary']] <- list_de_results[['Summary']]
list_de_results.significant.protein_coding[['Summary']] <- list_de_results.protein_coding[['Summary']]

for (i in 2:length(list_de_results)) {
  # All expressed
  df.de <- list_de_results[[i]]
  df.de$F <- NULL
  cnames <- colnames(df.de); cnames <- gsub('logFC', 'log2FC', cnames); colnames(df.de) <- cnames
  list_de_results[[i]] <- df.de
  df.de.sig <- df.de[df.de$FDR < FDR_threshold & abs(df.de$log2FC) > logFC_threshold, ]
  list_de_results.significant[[i]] <- df.de.sig

  # Coding only
  df.de <- list_de_results.protein_coding[[i]]
  df.de$F <- NULL
  cnames <- colnames(df.de); cnames <- gsub('logFC', 'log2FC', cnames); colnames(df.de) <- cnames
  list_de_results.protein_coding[[i]] <- df.de
  df.de.sig <- df.de[df.de$FDR < FDR_threshold & abs(df.de$log2FC) > logFC_threshold, ]
  list_de_results.significant.protein_coding[[i]] <- df.de.sig
}

## ---- Sheet naming and output (unchanged)
sheet_names = names(list_de_results)
sheet_names= gsub('Cntrl', 'CNT', sheet_names)
sheet_names= gsub('__vs__', '.vs.', sheet_names)
sheet_names= gsub('beta', 'B', sheet_names)

names(list_de_results)                         = sheet_names
names(list_de_results.protein_coding)          = sheet_names
names(list_de_results.significant)             = sheet_names
names(list_de_results.significant.protein_coding) = sheet_names

result_xlsx = paste0(dea_result_dir, '/DE_Results.All_expressed_coding_genes.xlsx')
write.xlsx(list_de_results.protein_coding, result_xlsx)

result_xlsx = paste0(dea_result_dir, '/DE_Results.All_expressed_total_genes.xlsx')
write.xlsx(list_de_results, result_xlsx)

result_xlsx = paste0(dea_result_dir, '/DE_Results.Significantly_DE_coding_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.xlsx')
write.xlsx(list_de_results.significant.protein_coding, result_xlsx)

result_xlsx = paste0(dea_result_dir, '/DE_Results.Significantly_DE_total_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.xlsx')
write.xlsx(list_de_results.significant, result_xlsx)

result_rds = paste0(dea_result_dir, '/DE_Results.All_expressed_coding_genes.rds')
saveRDS(list_de_results.protein_coding, result_rds)

result_rds = paste0(dea_result_dir, '/DE_Results.All_expressed_total_genes.rds')
saveRDS(list_de_results, result_rds)

result_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_coding_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.rds')
saveRDS(list_de_results.significant.protein_coding, result_rds)

result_rds = paste0(dea_result_dir, '/DE_Results.Significantly_DE_total_genes.FDR_',FDR_threshold, '.FC_',FC_threshold,'.rds')
saveRDS(list_de_results.significant, result_rds)
