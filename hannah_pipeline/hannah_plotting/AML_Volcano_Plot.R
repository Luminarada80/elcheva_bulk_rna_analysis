library(dplyr)
library(ggplot2)
library(ggrepel)
library(ragg)

# load
data <- readRDS("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/DIFFERENTIAL_EXPRESSION/2025_01/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest/DE_Results.All_expressed_coding_genes.rds")
data <- data$K562_shCNTL_vs_K562_shNYN3


#KG1A_shCNTL_vs_KG1A_shNYN4
#K562_shCNTL_vs_K562_shNYN3
#K562_shCNTL_vs_K562_shNYN4
#SKNAS_shCNTL_vs_SKNAS_shNYN4

# ensure unique column names, and DO NOT blindly rename column 3
data <- as.data.frame(data)
data <- data[, !duplicated(names(data))]

# If your log fold change column is named logFC, normalize it:
if (!"log2FC" %in% names(data) && "logFC" %in% names(data)) {
  data <- dplyr::rename(data, log2FC = logFC)
}

# If your significance column is called something else, point FDR to it.
# (Comment this if you already have 'FDR'.)
if (!"FDR" %in% names(data)) {
  # try common alternatives; pick the first match that exists
  cand <- c("FDR_log2", "neglog2FDR", "FDR")
  have <- cand[cand %in% names(data)]
  stopifnot(length(have) >= 1)
  data <- dplyr::rename(data, FDR = !!have[1])
}

# Optional: auto-convert if FDR looks like raw 0..1 values instead of -log2(FDR)
if (max(data$FDR, na.rm = TRUE) <= 1) {
  data$FDR <- -log2(pmax(data$FDR, .Machine$double.xmin))
}

# thresholds
th <- log2(1.5)

# DE flags + labels
data$diffexpressed <- "NO"
data$diffexpressed[data$log2FC >  th & data$FDR > 2] <- "UP"
data$diffexpressed[data$log2FC < -th & data$FDR > 2] <- "DOWN"
data$delabel <- ifelse(data$gene_symbol %in% c("CMYA5","H3C10","SCRN1","SARM1","OPRD1"),
                       data$gene_symbol, NA)

theme_set(theme_classic(base_size = 20) +
            theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), color = 'black'),
                  axis.title.x = element_text(face = "bold", margin = margin(20,0,0,0), color = 'black'),
                  plot.title   = element_text(hjust = 0.5)))

p <- ggplot(data, aes(x = log2FC, y = FDR, col = diffexpressed)) +
  geom_vline(xintercept = c(-th, th), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = 2,          col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  # draw black points for labeled genes without inheriting color mapping
  geom_point(data = subset(data, !is.na(delabel)),
             inherit.aes = FALSE,
             aes(x = log2FC, y = FDR),
             size = 2, col = "black") +
  # labels MUST provide label aesthetic
  geom_text_repel(data = subset(data, !is.na(delabel)),
                  aes(x = log2FC, y = FDR, label = delabel),
                  max.overlaps = Inf, color = "black") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "grey", "UP" = "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(color = 'Legend',
       x = expression("log"[2]*"FC"),
       y = expression("-log"[2]*"FDR")) +
  ggtitle('K562_shCNTL_vs_K562_shNYN3')

p

# save (no X11 required)
ggsave("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Sig_DE_Volcano_K562_shCNTL_vs_K562_shNYN3.png",
       p, device = ragg::agg_png, width = 9, height = 7, units = "in", dpi = 300, bg = "white")

