library(dplyr)
library(ggplot2)
library(ggrepel)
library(ragg)

# Load the list of DE results (one element per comparison)
data_list <- readRDS("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/DIFFERENTIAL_EXPRESSION/2025_01/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest/DE_Results.All_expressed_coding_genes.rds")

sample_names <- c(
  "KG1A_shCNTL_vs_KG1A_shNYN4",
  "K562_shCNTL_vs_K562_shNYN3",
  "K562_shCNTL_vs_K562_shNYN4",
  "SKNAS_shCNTL_vs_SKNAS_shNYN4"
)

# Global theme once
theme_set(
  theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), color = 'black'),
      axis.title.x = element_text(face = "bold", margin = margin(20,0,0,0), color = 'black'),
      plot.title   = element_text(hjust = 0.5)
    )
)

for (sample_name in sample_names) {
  # Check sample exists
  if (!sample_name %in% names(data_list)) {
    message("Skipping missing sample: ", sample_name)
    next
  }
  
  df <- data_list[[sample_name]]
  if (is.null(df) || nrow(df) == 0) {
    message("No rows for: ", sample_name)
    next
  }
  
  # Ensure unique column names (drop duplicates)
  df <- as.data.frame(df)
  df <- df[, !duplicated(names(df)), drop = FALSE]
  
  # Ensure log2FC column
  if (!"log2FC" %in% names(df)) {
    if ("logFC" %in% names(df)) {
      df <- dplyr::rename(df, log2FC = logFC)
    } else {
      message("Missing log2FC/logFC in: ", sample_name, " — skipping")
      next
    }
  }
  
  # Ensure FDR column (or a stand-in)
  if (!"FDR" %in% names(df)) {
    candidates <- c("FDR_log2", "neglog2FDR", "padj", "adj.P.Val", "qval", "FDR")
    have <- candidates[candidates %in% names(df)]
    if (length(have) == 0) {
      message("Missing FDR-like column in: ", sample_name, " — skipping")
      next
    }
    df <- dplyr::rename(df, FDR = !!have[1])
  }
  
  # Convert FDR: if looks like raw 0..1 p/q-values, convert to -log2
  if (is.numeric(df$FDR) && max(df$FDR, na.rm = TRUE) <= 1) {
    df$FDR <- -log2(pmax(df$FDR, .Machine$double.xmin))
  }
  
  # Gene symbol presence
  if (!"gene_symbol" %in% names(df)) {
    message("Missing gene_symbol in: ", sample_name, " — skipping labels (plot will still render)")
    df$gene_symbol <- NA_character_
  }
  
  # Thresholds
  th <- log2(1.5)
  yth <- -log10(0.05)
  
  # Flags & labels
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FC >  th & df$FDR > yth] <- "UP"
  df$diffexpressed[df$log2FC < -th & df$FDR > yth] <- "DOWN"
  
  # Example labels (adjust or remove as needed)
  df$delabel <- ifelse(df$gene_symbol %in% c("CMYA5","H3C10","SCRN1","SARM1","OPRD1"),
                       df$gene_symbol, NA)
  
  # Plot
  p <- ggplot(df, aes(x = log2FC, y = FDR, col = diffexpressed)) +
    geom_vline(xintercept = c(-th, th), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = yth,          col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    geom_point(
      data = subset(df, !is.na(delabel)),
      inherit.aes = FALSE,
      aes(x = log2FC, y = FDR),
      size = 2, col = "black"
    ) +
    geom_text_repel(
      data = subset(df, !is.na(delabel)),
      aes(x = log2FC, y = FDR, label = delabel),
      max.overlaps = Inf, color = "black"
    ) +
    scale_color_manual(
      values = c("DOWN" = "#00AFBB", "NO" = "grey", "UP" = "#bb0c00"),
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"),
         y = expression("-log"[10]*"FDR")) +
    ggtitle(sample_name)
  
  # Save (no X11 needed)
  out_png <- sprintf("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Sig_DE_Volcano_%s.png", sample_name)
  ggsave(out_png, p, device = ragg::agg_png, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
  
  message("Saved: ", out_png)
}
