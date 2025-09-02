library(dplyr)
library(pheatmap)

# Load the data (adjust the file path accordingly)
datafile <- "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/DIFFERENTIAL_EXPRESSION/2025_01/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest/DE_Results.Significantly_DE_coding_genes.FDR_0.05.FC_1.5.rds"
# datafile <- "/gpfs/Labs/Uzun/DATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/BULK_RNA_SEQ/GENE_EXPRESSION_MATRICES/2025_01/FPKM_Matrix.cleaned.no_annot.gene_symbols_as_rownames.rds"
# data <- readRDS("DE_Results.Significantly_DE_coding_genes.FDR_0.05.FC_1.5.rds")
data <- readRDS(datafile)
data <- data$SKNAS_shCNTL_vs_SKNAS_shNYN4# Adjust this if necessary
#K562_shCNTL_vs_K562_shNYN4
#SKNAS_shCNTL_vs_SKNAS_shNYN4


# Sort by FDR
sorted_data <- data[order(data$FDR), ]

# Extract upregulated and downregulated genes
upregulated <- subset(sorted_data, log2FC > 0)
downregulated <- subset(sorted_data, log2FC < 0)

# Combine top 10 upregulated and downregulated genes
top_10_upregulated <- head(upregulated, 10)
top_10_downregulated <- head(downregulated, 10)
combined_data <- rbind(top_10_upregulated, top_10_downregulated)

# Extract gene_symbol and the last 6 columns (expression values)
gene_symbols <- combined_data$gene_symbol
last_6_columns <- combined_data[, (ncol(combined_data) - 5):ncol(combined_data)]
new_data <- data.frame(gene_symbol = gene_symbols, last_6_columns)

# Convert data frame to numerical matrix (removing gene_symbol column)
matrix_data <- as.matrix(new_data[, -1])
rownames(matrix_data) <- new_data$gene_symbol

# Perform log transformation on the matrix (adding a small value to avoid log(0))
log_matrix <- log(matrix_data + 1)

# Plot heatmap
plot <- pheatmap(log_matrix, border_color = "black",
         color = colorRampPalette(c("white", "red"))(100), 
         fontsize = 20, number_color = "black",
         main = "FPKM: SKNAS_shCNTL_vs_SKNAS_shNYN4",
         breaks = seq(min(log_matrix, na.rm = TRUE), max(log_matrix, na.rm = TRUE), length.out = 100),
         treeheight_col = 0, treeheight_row = 0,
         cluster_cols = FALSE, cluster_rows = FALSE)
ggsave("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/FPKM_Heatmap_SKNAS_shCNTL_vs_SKNAS_shNYN4.png", plot, width=10, height=10, units="in")

