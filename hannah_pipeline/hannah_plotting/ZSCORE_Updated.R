library(dplyr)
library(pheatmap)

# Load the data (adjust the file path accordingly)
data <- readRDS("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/DIFFERENTIAL_EXPRESSION/2025_01/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest/DE_Results.Significantly_DE_total_genes.FDR_0.05.FC_1.5.Split.rds")
data <- data$K562_shCNTL_vs_K562_shNYN4.Down # Adjust this if necessary

# Sort by FDR
sorted_data <- data[order(data$FDR), ]

# Remove mitochondrial genes (genes starting with "MT-" or "mt-")
non_mitochondrial_data <- sorted_data[!grepl("^MT-|^mt-", sorted_data$gene_symbol), ]

# Extract upregulated and downregulated genes
upregulated <- subset(non_mitochondrial_data, log2FC > 0)
downregulated <- subset(non_mitochondrial_data, log2FC < 0)

# Combine top 10 upregulated and downregulated geneshttp://127.0.0.1:18011/graphics/dfe15331-4770-49a5-b997-7ee1e58b0f42.png
top_10_upregulated <- head(upregulated, 10)
top_10_downregulated <- head(downregulated, 10)
combined_data <- rbind(top_10_upregulated, top_10_downregulated)

# Extract gene_symbol and the last 6 columns (expression values)
gene_symbols <- combined_data$gene_symbol
last_6_columns <- combined_data[, (ncol(combined_data) - 6):ncol(combined_data)]
new_data <- data.frame(gene_symbol = gene_symbols, last_6_columns)

# Convert data frame to numerical matrix (removing gene_symbol column)
matrix_data <- as.matrix(new_data[, -1])
rownames(matrix_data) <- new_data$gene_symbol

# Calculate z-score (row-wise)
z_score_matrix <- t(scale(t(matrix_data)))

# Set breaks for z-score at -2 and 2
z_breaks <- seq(-2, 2, length.out = 101)

# Plot heatmap
pheatmap(z_score_matrix, 
         border_color = "black",
         color = colorRampPalette(c("blue", "white"))(100), #, "red"
         fontsize = 20, number_color = "black",
         main = "BM: Z-Score: K562_shCNTL_vs_K562_shNYN4 Down",
         breaks = z_breaks,
         treeheight_col = 0, treeheight_row = 0,
         cluster_cols = FALSE, cluster_rows = FALSE)
ggsave("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Z_Score_Heatmap_K562_shCNTL_vs_K562_shNYN4.png", plot, width=10, height=10, units="in")

