# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ragg)

# ---- READ YOUR DATA ----
# If you're pasting from a file, use read.table("path", header=FALSE, fill=TRUE)
# Here I'll assume you already have a data.frame 'df' with 4 cols like the paste below:
df <- read.table(text='
K562_shCNTL_Rep1_R1.fastq.gz 517052 K562_shCNTL_Rep1_R1_paired.fastq.gz 258849
K562_shCNTL_Rep1_R2.fastq.gz 481704 K562_shCNTL_Rep1_R1_unpaired.fastq.gz 120855
K562_shCNTL_Rep2_R1.fastq.gz 484510 K562_shCNTL_Rep1_R2_paired.fastq.gz 262105
K562_shCNTL_Rep2_R2.fastq.gz 446212 K562_shCNTL_Rep1_R2_unpaired.fastq.gz 990
K562_shCNTL_Rep3_R1.fastq.gz 470147 K562_shCNTL_Rep2_R1_paired.fastq.gz 249878
K562_shCNTL_Rep3_R2.fastq.gz 445286 K562_shCNTL_Rep2_R1_unpaired.fastq.gz 118586
K562_shNYN3_Rep1_R1.fastq.gz 493775 K562_shCNTL_Rep2_R2_paired.fastq.gz 246131
K562_shNYN3_Rep1_R2.fastq.gz 468233 K562_shCNTL_Rep2_R2_unpaired.fastq.gz 979
K562_shNYN3_Rep2_R1.fastq.gz 418593 K562_shCNTL_Rep3_R1_paired.fastq.gz 241359
K562_shNYN3_Rep2_R2.fastq.gz 382159 K562_shCNTL_Rep3_R1_unpaired.fastq.gz 110403
K562_shNYN3_Rep3_R1.fastq.gz 471766 K562_shCNTL_Rep3_R2_paired.fastq.gz 240556
K562_shNYN3_Rep3_R2.fastq.gz 420865 K562_shCNTL_Rep3_R2_unpaired.fastq.gz 950
K562_shNYN4_Rep1_R1.fastq.gz 477315 K562_shNYN3_Rep1_R1_paired.fastq.gz 248046
K562_shNYN4_Rep1_R2.fastq.gz 451345 K562_shNYN3_Rep1_R1_unpaired.fastq.gz 113845
K562_shNYN4_Rep2_R1.fastq.gz 516728 K562_shNYN3_Rep1_R2_paired.fastq.gz 253743
K562_shNYN4_Rep2_R2.fastq.gz 470662 K562_shNYN3_Rep1_R2_unpaired.fastq.gz 945
K562_shNYN4_Rep3_R1.fastq.gz 417553 K562_shNYN3_Rep2_R1_paired.fastq.gz 209953
K562_shNYN4_Rep3_R2.fastq.gz 394961 K562_shNYN3_Rep2_R1_unpaired.fastq.gz 100857
KG1A_shCNTL_Rep1_S18_R1.fastq.gz 452434 K562_shNYN3_Rep2_R2_paired.fastq.gz 209169
KG1A_shCNTL_Rep1_S18_R2.fastq.gz 408898 K562_shNYN3_Rep2_R2_unpaired.fastq.gz 837
KG1A_shCNTL_Rep2_S20_R1.fastq.gz 475089 K562_shNYN3_Rep3_R1_paired.fastq.gz 242575
KG1A_shCNTL_Rep2_S20_R2.fastq.gz 456464 K562_shNYN3_Rep3_R1_unpaired.fastq.gz 107390
KG1A_shCNTL_Rep3_S22_R1.fastq.gz 426666 K562_shNYN3_Rep3_R2_paired.fastq.gz 237541
KG1A_shCNTL_Rep3_S22_R2.fastq.gz 412221 K562_shNYN3_Rep3_R2_unpaired.fastq.gz 898
KG1A_shNYN4_Rep1_S19_R1.fastq.gz 363565 K562_shNYN4_Rep1_R1_paired.fastq.gz 245542
KG1A_shNYN4_Rep1_S19_R2.fastq.gz 329750 K562_shNYN4_Rep1_R1_unpaired.fastq.gz 112769
KG1A_shNYN4_Rep2_S21_R1.fastq.gz 432177 K562_shNYN4_Rep1_R2_paired.fastq.gz 245323
KG1A_shNYN4_Rep2_S21_R2.fastq.gz 431807 K562_shNYN4_Rep1_R2_unpaired.fastq.gz 1002
KG1A_shNYN4_Rep3_S23_R1.fastq.gz 556052 K562_shNYN4_Rep2_R1_paired.fastq.gz 264785
KG1A_shNYN4_Rep3_S23_R2.fastq.gz 517568 K562_shNYN4_Rep2_R1_unpaired.fastq.gz 120927
SKNAS_shCNTL_Rep1_R1.fastq.gz 373573 K562_shNYN4_Rep2_R2_paired.fastq.gz 261129
SKNAS_shCNTL_Rep1_R2.fastq.gz 355798 K562_shNYN4_Rep2_R2_unpaired.fastq.gz 996
SKNAS_shCNTL_Rep2_R1.fastq.gz 373231 K562_shNYN4_Rep3_R1_paired.fastq.gz 212972
SKNAS_shCNTL_Rep2_R2.fastq.gz 347718 K562_shNYN4_Rep3_R1_unpaired.fastq.gz 104190
SKNAS_shCNTL_Rep3_R1.fastq.gz 480186 K562_shNYN4_Rep3_R2_paired.fastq.gz 209898
SKNAS_shCNTL_Rep3_R2.fastq.gz 442038 K562_shNYN4_Rep3_R2_unpaired.fastq.gz 736
SKNAS_shNYN4_Rep1_R1.fastq.gz 422452 KG1A_shCNTL_Rep1_S18_R1_paired.fastq.gz 229155
SKNAS_shNYN4_Rep1_R2.fastq.gz 369729 KG1A_shCNTL_Rep1_S18_R1_unpaired.fastq.gz 104033
SKNAS_shNYN4_Rep2_R1.fastq.gz 420201 KG1A_shCNTL_Rep1_S18_R2_paired.fastq.gz 228258
SKNAS_shNYN4_Rep2_R2.fastq.gz 408932 KG1A_shCNTL_Rep1_S18_R2_unpaired.fastq.gz 938
SKNAS_shNYN4_Rep3_R1.fastq.gz 448010 KG1A_shCNTL_Rep2_S20_R1_paired.fastq.gz 238917
SKNAS_shNYN4_Rep3_R2.fastq.gz 437231 KG1A_shCNTL_Rep2_S20_R1_unpaired.fastq.gz 109880
       KG1A_shCNTL_Rep2_S20_R2_paired.fastq.gz 244047
       KG1A_shCNTL_Rep2_S20_R2_unpaired.fastq.gz 996
       KG1A_shCNTL_Rep3_S22_R1_paired.fastq.gz 212131
       KG1A_shCNTL_Rep3_S22_R1_unpaired.fastq.gz 100408
       KG1A_shCNTL_Rep3_S22_R2_paired.fastq.gz 215895
       KG1A_shCNTL_Rep3_S22_R2_unpaired.fastq.gz 832
       KG1A_shNYN4_Rep1_S19_R1_paired.fastq.gz 187111
       KG1A_shNYN4_Rep1_S19_R1_unpaired.fastq.gz 84308
       KG1A_shNYN4_Rep1_S19_R2_paired.fastq.gz 183688
       KG1A_shNYN4_Rep1_S19_R2_unpaired.fastq.gz 739
       KG1A_shNYN4_Rep2_S21_R1_paired.fastq.gz 242889
       KG1A_shNYN4_Rep2_S21_R1_unpaired.fastq.gz 108907
       KG1A_shNYN4_Rep2_S21_R2_paired.fastq.gz 238841
       KG1A_shNYN4_Rep2_S21_R2_unpaired.fastq.gz 865
       KG1A_shNYN4_Rep3_S23_R1_paired.fastq.gz 278750
       KG1A_shNYN4_Rep3_S23_R1_unpaired.fastq.gz 129355
       KG1A_shNYN4_Rep3_S23_R2_paired.fastq.gz 284296
       KG1A_shNYN4_Rep3_S23_R2_unpaired.fastq.gz 1041
       SKNAS_shCNTL_Rep1_R1_paired.fastq.gz 189253
       SKNAS_shCNTL_Rep1_R1_unpaired.fastq.gz 84399
       SKNAS_shCNTL_Rep1_R2_paired.fastq.gz 188560
       SKNAS_shCNTL_Rep1_R2_unpaired.fastq.gz 837
       SKNAS_shCNTL_Rep2_R1_paired.fastq.gz 191175
       SKNAS_shCNTL_Rep2_R1_unpaired.fastq.gz 80842
       SKNAS_shCNTL_Rep2_R2_paired.fastq.gz 188519
       SKNAS_shCNTL_Rep2_R2_unpaired.fastq.gz 783
       SKNAS_shCNTL_Rep3_R1_paired.fastq.gz 235965
       SKNAS_shCNTL_Rep3_R1_unpaired.fastq.gz 112044
       SKNAS_shCNTL_Rep3_R2_paired.fastq.gz 241344
       SKNAS_shCNTL_Rep3_R2_unpaired.fastq.gz 856
       SKNAS_shNYN4_Rep1_R1_paired.fastq.gz 209939
       SKNAS_shNYN4_Rep1_R1_unpaired.fastq.gz 94577
       SKNAS_shNYN4_Rep1_R2_paired.fastq.gz 212960
       SKNAS_shNYN4_Rep1_R2_unpaired.fastq.gz 873
       SKNAS_shNYN4_Rep2_R1_paired.fastq.gz 206683
       SKNAS_shNYN4_Rep2_R1_unpaired.fastq.gz 101643
       SKNAS_shNYN4_Rep2_R2_paired.fastq.gz 213203
       SKNAS_shNYN4_Rep2_R2_unpaired.fastq.gz 789
       SKNAS_shNYN4_Rep3_R1_paired.fastq.gz 219876
       SKNAS_shNYN4_Rep3_R1_unpaired.fastq.gz 111915
       SKNAS_shNYN4_Rep3_R2_paired.fastq.gz 229704
       SKNAS_shNYN4_Rep3_R2_unpaired.fastq.gz 784
', header = FALSE, fill = TRUE, col.names = paste0("V", 1:4))

# ---- CLEAN & COMPUTE ----
# Fill down missing V1/V2 from previous non-NA
df <- df %>%
  mutate(across(c(V1, V2), ~na_if(., ""))) %>%
  tidyr::fill(V1, V2, .direction = "down")

# Keep only "paired" rows in V3 (these are what you want in the numerator)
paired <- df %>% 
  filter(!is.na(V3), str_detect(V3, "_paired\\.fastq\\.gz$"))

# Compute percentage per row: paired_after / raw
paired <- paired %>%
  mutate(V2 = as.numeric(V2),
         V4 = as.numeric(V4),
         Percentage = 100 * V4 / V2)

# (Optional) If you want a single bar per *sample base* (collapse R1/R2):
# Extract a base sample ID without _R1/_R2
paired <- paired %>%
  mutate(Sample = str_replace(V1, "_R[12]\\.fastq\\.gz$", ""))

# Aggregate to one percentage per Sample by averaging R1 and R2
by_sample <- paired %>%
  group_by(Sample) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

# Mean across samples
mean_row <- tibble(Sample = "Mean", Percentage = round(mean(by_sample$Percentage, na.rm = TRUE), 2))

plot_df <- bind_rows(by_sample, mean_row)

# Order samples (optional: alphabetical, or keep original)
plot_df <- plot_df %>%
  mutate(Sample = factor(Sample, levels = c(sort(setdiff(Sample, "Mean")), "Mean")))

# ---- PLOT ----
p <- ggplot(plot_df, aes(x = Sample, y = Percentage)) +
  geom_col(fill = "skyblue") +
  geom_text(aes(label = sprintf("%.0f", Percentage)),
            vjust = 1.4, size = 3, angle = 0, hjust=0.5) +
  labs(title = "Percent of FASTQ Reads Retained After Trimming (paired/total)",
       x = "Sample", y = "Percent (%)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "grey95", colour = NA))

# Print to screen (if interactive)
p

# ---- SAVE HEADLESS-SAFE ----
ggsave(
  filename = "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Trimming_Percent_Reads_Retained_Boxplot.png",
  plot = p, device = ragg::agg_png, width = 8, height = 5, units = "in", dpi = 300, bg = "white"
)

