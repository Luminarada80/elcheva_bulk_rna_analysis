library(ggplot2)
library(dplyr)
library(tidyr)

data <- readRDS("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/DIFFERENTIAL_EXPRESSION/2025_01/FPKM_GT_ONE_IN_AT_LEAST_TWO_SAMPLES_IN_EITHER_GROUP/glmQLFTest/DE_Results.Significantly_DE_total_genes.FDR_0.05.FC_1.5.Split.rds")

data <- data[-1]

data <- lapply(data, function(df) subset(df, PValue <= 0.01))

data <- lapply(data, function(df) subset(df, abs(log2FC) >= 1.5))

data <- lapply(data, function(df) subset(df, G1_MEAN_FPKM >= 5))

data <- lapply(data, function(df) subset(df, G2_MEAN_FPKM >= 5))

upregulated <- lapply(data, function(df) subset(df, log2FC > 0))

downregulated <- lapply(data, function(df) subset(df, log2FC < 0))

upregulated_count <- lapply(upregulated, function(df) nrow(df))

downregulated_count <- lapply(downregulated, function(df) nrow(df))

combined_counts <- data.frame(
  Comparison = names(upregulated_count),
  Upregulated_Count = unlist(upregulated_count),
  Downregulated_Count = unlist(downregulated_count)
)

# Create your data frame
combined_counts <- data.frame(
  Comparison = c("KG1A_shCTRL_vs_KG1A_shNYN4", "KG1A_shCTRL_vs_KG1A_shNYN4",
                 "K562_shCTRL_vs_K562_shNYN3", "K562_shCTRL_vs_K562_shNYN3",
                 "K562_shCTRL_vs_K562_shNYN4", "K562_shCTRL_vs_K562_shNYN4",
                 "SKNAS_shCTRL_vs_SKNAS_shNYN4", "SKNAS_shCTRL_vs_SKNAS_shNYN4"),
  Type = c("Upregulated", "Downregulated", 
           "Upregulated", "Downregulated", 
           "Upregulated", "Downregulated",
           "Upregulated", "Downregulated"),
  Count = c(0, 0, 0, 0, 0, 3, 0, 9)
)

# Plot using ggplot2
plot <- ggplot(combined_counts, aes(x=Comparison, y=Count, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("Upregulated"="red", "Downregulated"="blue")) +
  geom_text(aes(label=Count), 
            position=position_dodge(width=0.9), 
            vjust=-0.5, size=3) +
  theme_minimal() +
  labs(title="Upregulated and Downregulated Counts",
       x="Comparison",
       y="Count",
       fill="Regulation Type") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 65, hjust = 1))

plot

ggsave(
  filename = "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Up_Down_Reg_Sig_DE_Genes_by_Comparison.png",
  plot = plot, device = ragg::agg_png, width = 8, height = 5, units = "in", dpi = 300, bg = "white"
)
