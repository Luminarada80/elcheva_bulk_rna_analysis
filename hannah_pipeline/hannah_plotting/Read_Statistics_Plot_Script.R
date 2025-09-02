
library(ggplot2)

# Read the data
read_stats <- read.table("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/READ_STATS/MERGED/Read_Statistics.2025_01.txt", header = TRUE)
read_stats <- read_stats[-which(names(read_stats) == "Unique")]
read_stats$Total <- rowSums(read_stats[, -1])

# Calculate percentages
read_stats$Intergenic_Percent <- (read_stats$Intergenic / read_stats$Total) * 100
read_stats$Intron_Percent <- (read_stats$Intron / read_stats$Total) * 100
read_stats$Junction_Percent <- (read_stats$Junction / read_stats$Total) * 100
read_stats$Exon_Percent <- (read_stats$Exon / read_stats$Total) * 100

# Subset the data
read_stats <- subset(read_stats, select = -c(Exon, Intron, Intergenic, Junction, Total))

# Reshape the data into long format
read_stats_long <- reshape(read_stats,
                           varying = list(names(read_stats)[-1]),
                           v.names = "Count",
                           timevar = "Category",
                           times = names(read_stats)[-1],
                           idvar = "Sample",
                           direction = "long")
read_stats_long$Category <- factor(read_stats_long$Category, levels = c("Intergenic_Percent", "Intron_Percent", "Junction_Percent", "Exon_Percent"))

# Create the plot
plot <- ggplot(read_stats_long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Read Statistics by Sample",
       x = "Sample",
       y = "Read Count %",
       fill = "Category") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),

        axis.text.x = element_text(angle = 65, hjust = 1))

# Calculate means
means <- colMeans(read_stats[, -1], na.rm = TRUE)

# Create data frame for means
means_df <- data.frame(Sample = rep("Mean", length(means)), Category = names(means), Count = means)
means_df$Category <- factor(means_df$Category,
                            levels = c("Intergenic_Percent","Intron_Percent","Junction_Percent","Exon_Percent"))

# Add means to the plot
plot <- plot + geom_bar(data = means_df, aes(x = Sample, y = Count, fill = Category), stat = "identity", position = "stack")

# Show the plot
plot

ggsave("/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/PLOTS/Read_Statistics_Stacked_Barplot.png", plot, width=8, height=6, units="in")
