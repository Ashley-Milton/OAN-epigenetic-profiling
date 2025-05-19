#Set working directory
setwd("/path/to/working/dir")

#Load packages
library(BiocManager)
library(GenomicRanges)
library(InteractionSet)
library(tidyverse)
library(reshape2)
library(ggthemes)
library(viridis)
library(readr)
library(data.table)

#Define binsize variable
binsize <- "500kb"
bin_size <- as.numeric(gsub("kb", "000", binsize))

#Load male and female plotting data (created in previous sex-specific inter-chromosomal interaction plot scripts)
data_F <- fread("data/longDatafemale_500kb_zerowindows.txt")
data_M <- fread("data/longDatamale_500kb_zerowindows.txt")

#Combine data by Var1 and Var2
data_combined <- merge(data_F, data_M, by = c("Var1", "Var2"), suffixes = c("_F", "_M"))

#Add ratio column
data_combined$ratio <- data_combined$value_F / data_combined$value_M

#Filter data_combined to remove rows where Var1 is equal to Var2
filtered_data_combined <- filter(data_combined, Var1 != Var2)

#Calculate lim1 and lim2
min_value <- min(filtered_data_combined$ratio)
max_value <- max(filtered_data_combined$ratio)

#Auto limits
lim1 <- floor(min_value / 0.0005) * 0.0005
lim2 <- ceiling(max_value / 0.0005) * 0.0005

#Overriding with manual limits
lim1 <- 0.4
lim2 <- 1.6

midp <- lim1 + ((lim2 - lim1) / 2)

#Plot data
plot <- ggplot(data_combined, aes(Var2, Var1, fill = ratio)) + 
  theme_few() +
  ggtitle(paste("Inter-chromosomal ratio for platypus (", binsize, " bins)", sep = "")) +
  geom_tile(colour = "grey", size = 0.5) +
  scale_fill_gradient2(name = expression("Ratio\nInteractions"), low = "navyblue", high = "#FF4500", mid = "white", midpoint = midp, limits = c(lim1, lim2)) +
  scale_x_discrete(name = "Chromosomes", limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X1", "X2", "X3", "X4", "X5")) + 
  scale_y_discrete(name = "Chromosomes", limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X1", "X2", "X3", "X4", "X5"))

plot

#Save plot
ggsave(paste0("plots/Platypus_", binsize, "_lims_", lim1, "_", lim2, "_zerowindows_differential.pdf"), plot = plot, width = 10, height = 8, units = "in", dpi = 800)

