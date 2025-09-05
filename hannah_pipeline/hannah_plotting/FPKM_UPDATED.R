library(dplyr)
library(pheatmap)

datafile <- "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/DIFFERENTIAL_EXPRESSION/2025_01/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest/DE_Results.Significantly_DE_coding_genes.FDR_0.05.FC_1.5.rds"
data_list <- readRDS(datafile)  # keep the full list separate

sample_names <- c(
  "KG1A_shCNTL_vs_KG1A_shNYN4",
  "K562_shCNTL_vs_K562_shNYN3",
  "K562_shCNTL_vs_K562_shNYN4",
  "SKNAS_shCNTL_vs_SKNAS_shNYN4"
)

for (sample_name in sample_names) {
  # skip if the sample is missing
  if (!sample_name %in% names(data_list)) {
    message("Skipping missing sample: ", sample_name)
    next
  }
  
  df <- data_list[[sample_name]]
  if (is.null(df) || nrow(df) == 0) {
    message("No rows for: ", sample_name)
    next
  }
  
  # Sort by FDR
  sorted_df <- df[order(df$FDR), ]
  
  # Split by direction
  upregulated <- subset(sorted_df, log2FC > 0)
  downregulated <- subset(sorted_df, log2FC < 0)
  
  # Top 10 from each
  top_10_up <- head(upregulated, 10)
  top_10_down <- head(downregulated, 10)
  combined <- rbind(top_10_up, top_10_down)
  
  if (nrow(combined) == 0) {
    message("No DE genes to plot for: ", sample_name)
    next
  }
  
  # Take the last 6 columns as expression (adjust if needed)
  expr <- combined[, (ncol(combined) - 5):ncol(combined)]
  
  # Ensure numeric matrix
  expr <- as.data.frame(lapply(expr, function(x) as.numeric(as.character(x))))
  mat <- as.matrix(expr)
  rownames(mat) <- combined$gene_symbol
  
  # Log transform
  log_mat <- log(mat + 1)
  
  # Either let pheatmap pick breaks OR set 101 breaks for 100 colors
  cols <- colorRampPalette(c("blue", "white", "red"))(100)
  brks <- seq(min(log_mat, na.rm = TRUE), max(log_mat, na.rm = TRUE), length.out = 101)
  
  out_png <- paste0("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/FPKM_Heatmap_", sample_name, ".png")
  
  # Save directly via pheatmap
  pheatmap(
    log_mat,
    color = cols,
    breaks = brks,
    border_color = "black",
    fontsize = 20,
    main = paste0("FPKM: ", sample_name),
    treeheight_col = 0,
    treeheight_row = 0,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    filename = out_png,      # <- saves to file
    width = 10, height = 10  # inches
  )
}
