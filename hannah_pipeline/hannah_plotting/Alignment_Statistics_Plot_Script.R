
library(ggplot2)

# Read the data
alignment_data <- read.table("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/STAR_ALIGNMENT_STATISTICS/STAR_Alignment_Statistics.2025_01.txt", header = TRUE, sep = "\t")

# Remove the % sign and convert to numeric
alignment_data$Unique_Alignment_Rate <- as.numeric(gsub("%", "", alignment_data$Unique_Alignment_Rate))

# Calculate the mean of Unique_Alignment_Rate
mean_alignment_rate <- mean(alignment_data$Unique_Alignment_Rate)

# Round the mean to integer
mean_alignment_rate_rounded <- round(mean_alignment_rate)

# Create a new row for Mean_Unique_Alignment_Rate
mean_alignment_row <- data.frame(Sample = "Mean_Unique_Alignment_Rate", Unique_Alignment_Rate = mean_alignment_rate_rounded)

# Add Unique_Alignment_Rate_Integer column
alignment_data$Unique_Alignment_Rate_Integer <- round(alignment_data$Unique_Alignment_Rate)

# Create the plot
plot <- ggplot() +
  geom_bar(data = alignment_data, aes(x = Sample, y = Unique_Alignment_Rate), fill = "skyblue", stat = "identity") +
  geom_bar(data = mean_alignment_row, aes(x = Sample, y = Unique_Alignment_Rate), fill = "gold2", stat = "identity") +
  geom_text(data = alignment_data, aes(x = Sample, y = Unique_Alignment_Rate, label = paste0(Unique_Alignment_Rate_Integer)), vjust = -0.5, size = 4) +
  labs(title = "Alignment Rates for Different Samples",
       x = "Sample", y = "Alignment Rate") +
  theme(axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) +
  ylim(0, 100) +
  geom_text(data = mean_alignment_row, aes(x = Sample, y = Unique_Alignment_Rate, label = paste0(Unique_Alignment_Rate)), vjust = -0.5, size = 4) +
  theme(legend.position = "none")

plot

ggsave("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Alignment_Statistics_Plot.png", plot, width=8, height=6, units="in")

