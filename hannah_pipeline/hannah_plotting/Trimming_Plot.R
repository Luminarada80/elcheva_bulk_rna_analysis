# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ragg)

# ---- READ YOUR DATA ----
df <- read.table(text='
K562_shCNTL_Rep1_R1.fastq.gz	517052	K562_shCNTL_Rep1_R1_paired.fastq.gz	260296
K562_shCNTL_Rep1_R2.fastq.gz	481704	K562_shCNTL_Rep1_R1_unpaired.fastq.gz	354850
K562_shCNTL_Rep2_R1.fastq.gz	484510	K562_shCNTL_Rep1_R2_unpaired.fastq.gz	867
K562_shCNTL_Rep2_R2.fastq.gz	446212	K562_shCNTL_Rep2_R1_paired.fastq.gz	249514
K562_shCNTL_Rep3_R1.fastq.gz	470147	K562_shCNTL_Rep2_R1_unpaired.fastq.gz	344783
K562_shCNTL_Rep3_R2.fastq.gz	445286	K562_shCNTL_Rep3_R1_paired.fastq.gz	241594
K562_shNYN3_Rep1_R1.fastq.gz	493775	K562_shCNTL_Rep3_R1_unpaired.fastq.gz	327781
K562_shNYN3_Rep1_R2.fastq.gz	468233	K562_shCNTL_Rep3_R2_unpaired.fastq.gz	836
K562_shNYN3_Rep2_R1.fastq.gz	418593	K562_shNYN3_Rep1_R1_paired.fastq.gz	247615
K562_shNYN3_Rep2_R2.fastq.gz	382159	K562_shNYN3_Rep1_R1_unpaired.fastq.gz	334451
K562_shNYN3_Rep3_R1.fastq.gz	471766	K562_shNYN3_Rep1_R2_unpaired.fastq.gz	875
K562_shNYN3_Rep3_R2.fastq.gz	420865	K562_shNYN3_Rep2_R1_paired.fastq.gz	209404
K562_shNYN4_Rep1_R1.fastq.gz	477315	K562_shNYN3_Rep2_R1_unpaired.fastq.gz	289221
K562_shNYN4_Rep1_R2.fastq.gz	451345	K562_shNYN3_Rep3_R1_paired.fastq.gz	242916
K562_shNYN4_Rep2_R1.fastq.gz	516728	K562_shNYN3_Rep3_R1_unpaired.fastq.gz	325526
K562_shNYN4_Rep2_R2.fastq.gz	470662	K562_shNYN4_Rep1_R1_paired.fastq.gz	245628
K562_shNYN4_Rep3_R1.fastq.gz	417553	K562_shNYN4_Rep1_R1_unpaired.fastq.gz	334621
K562_shNYN4_Rep3_R2.fastq.gz	394961	K562_shNYN4_Rep2_R1_paired.fastq.gz	264640
KG1A_shCNTL_Rep1_S18_R1.fastq.gz	452434	K562_shNYN4_Rep2_R1_unpaired.fastq.gz	361603
KG1A_shCNTL_Rep1_S18_R2.fastq.gz	408898	K562_shNYN4_Rep3_R1_paired.fastq.gz	213026
KG1A_shCNTL_Rep2_S20_R1.fastq.gz	475089	K562_shNYN4_Rep3_R1_unpaired.fastq.gz	292635
KG1A_shCNTL_Rep2_S20_R2.fastq.gz	456464	KG1A_shCNTL_Rep1_S18_R1_paired.fastq.gz	229055
KG1A_shCNTL_Rep3_S22_R1.fastq.gz	426666	KG1A_shCNTL_Rep1_S18_R1_unpaired.fastq.gz	308471
KG1A_shCNTL_Rep3_S22_R2.fastq.gz	412221	KG1A_shCNTL_Rep2_S20_R1_paired.fastq.gz	238895
KG1A_shNYN4_Rep1_S19_R1.fastq.gz	363565	KG1A_shCNTL_Rep2_S20_R1_unpaired.fastq.gz	322134
KG1A_shNYN4_Rep1_S19_R2.fastq.gz	329750	KG1A_shCNTL_Rep3_S22_R1_paired.fastq.gz	212139
KG1A_shNYN4_Rep2_S21_R1.fastq.gz	432177	KG1A_shCNTL_Rep3_S22_R1_unpaired.fastq.gz	286348
KG1A_shNYN4_Rep2_S21_R2.fastq.gz	431807	KG1A_shNYN4_Rep1_S19_R1_paired.fastq.gz	187109
KG1A_shNYN4_Rep3_S23_R1.fastq.gz	556052	KG1A_shNYN4_Rep1_S19_R1_unpaired.fastq.gz	250937
KG1A_shNYN4_Rep3_S23_R2.fastq.gz	517568	KG1A_shNYN4_Rep2_S21_R1_paired.fastq.gz	242332
SKNAS_shCNTL_Rep1_R1.fastq.gz	373573	KG1A_shNYN4_Rep2_S21_R1_unpaired.fastq.gz	329858
SKNAS_shCNTL_Rep1_R2.fastq.gz	355798	KG1A_shNYN4_Rep3_S23_R1_paired.fastq.gz	279217
SKNAS_shCNTL_Rep2_R1.fastq.gz	373231	KG1A_shNYN4_Rep3_S23_R1_unpaired.fastq.gz	380491
SKNAS_shCNTL_Rep2_R2.fastq.gz	347718	SKNAS_shCNTL_Rep1_R1_paired.fastq.gz	189446
SKNAS_shCNTL_Rep3_R1.fastq.gz	480186	SKNAS_shCNTL_Rep1_R1_unpaired.fastq.gz	248180
SKNAS_shCNTL_Rep3_R2.fastq.gz	442038	SKNAS_shCNTL_Rep2_R1_paired.fastq.gz	191339
SKNAS_shNYN4_Rep1_R1.fastq.gz	422452	SKNAS_shCNTL_Rep2_R1_unpaired.fastq.gz	250179
SKNAS_shNYN4_Rep1_R2.fastq.gz	369729	SKNAS_shCNTL_Rep3_R1_paired.fastq.gz	236548
SKNAS_shNYN4_Rep2_R1.fastq.gz	420201	SKNAS_shCNTL_Rep3_R1_unpaired.fastq.gz	319108
SKNAS_shNYN4_Rep2_R2.fastq.gz	408932	SKNAS_shNYN4_Rep1_R1_paired.fastq.gz	209922
SKNAS_shNYN4_Rep3_R1.fastq.gz	448010	SKNAS_shNYN4_Rep1_R1_unpaired.fastq.gz	279281
SKNAS_shNYN4_Rep3_R2.fastq.gz	437231	SKNAS_shNYN4_Rep2_R1_paired.fastq.gz	207230
	SKNAS_shNYN4_Rep2_R1_unpaired.fastq.gz	282309
	SKNAS_shNYN4_Rep3_R1_paired.fastq.gz	219777
	SKNAS_shNYN4_Rep3_R1_unpaired.fastq.gz	305540

', header = FALSE, fill = TRUE, col.names = paste0("V", 1:4))

# ---- CLEAN & COMPUTE ----
# Treat blanks as NA and fill V1/V2 down
df <- df %>%
  mutate(
    V1 = na_if(V1, ""),
    V3 = na_if(V3, "")
  ) %>%
  fill(V1, V2, .direction = "down")


# Keep only the "paired" output rows in V3
paired <- df %>%
  filter(!is.na(V3) & str_detect(V3, "_paired\\.fastq\\.gz$"))

# Defensive numeric conversion
paired <- paired %>%
  mutate(
    V2 = suppressWarnings(as.numeric(V2)),
    V4 = suppressWarnings(as.numeric(V4))
  )

# Guard against divide-by-zero and NA
paired <- paired %>%
  mutate(
    Percentage = dplyr::if_else(!is.na(V2) & V2 > 0, 100 * V4 / V2, NA_real_)
  )

# Extract base sample ID (drop _R1/_R2.fastq.gz at end of V1)
paired <- paired %>%
  mutate(Sample = str_replace(V1, "_R[12]\\.fastq\\.gz$", ""))

# Aggregate to one value per Sample (average R1/R2)
by_sample <- paired %>%
  group_by(Sample) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

# Mean across samples
mean_val <- round(mean(by_sample$Percentage, na.rm = TRUE), 2)
plot_df <- bind_rows(by_sample, tibble(Sample = "Mean", Percentage = mean_val))

# Order: alphabetical samples, then "Mean"
plot_df <- plot_df %>%
  mutate(Sample = factor(Sample, levels = c(sort(setdiff(Sample, "Mean")), "Mean")))

# ---- PLOT ----
p <- ggplot(plot_df, aes(x = Sample, y = Percentage)) +
  geom_col(fill = "skyblue") +
  geom_hline(yintercept = mean_val, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.0f", Percentage)),
            vjust = 1.4, size = 3) +
  labs(
    title = "Percent of FASTQ Reads Retained After Trimming (paired / total)",
    x = "Sample", y = "Percent (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.background = element_rect(fill = "grey95", colour = NA)
  )

p

# ---- SAVE HEADLESS-SAFE ----
ggsave(
  filename = "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Trimming_Percent_Reads_Retained_Bar.png",
  plot = p, device = ragg::agg_png, width = 9, height = 5.5, units = "in",
  dpi = 300, bg = "white"
)
