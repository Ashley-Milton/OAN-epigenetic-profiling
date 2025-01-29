#Set working directory
setwd("/path/to/working/dir")

#Load libraries
library(BiocManager)
library(GenomicRanges)
library(InteractionSet)
library(tidyverse)
library(reshape2)
library(ggthemes)
library(viridis)
library(readr)

#Define binsize variable
binsize <- "500kb"
bin_size <- as.numeric(gsub("kb", "000", binsize))

#Read data for male
hic <- read.delim(paste0("data/Male_matrix_subset_2samplesonly_13_", binsize, "_bin_normalised_corrected_-2_3.tsv"), header = FALSE)

#Define the replacement rules
replacement_rules <- setNames(c(22:26, 27:31), c(paste0("X", 1:5), paste0("Y", 1:5)))

#Apply the replacement rules to columns 1 and 4
hic$V1 <- ifelse(hic$V1 %in% names(replacement_rules), replacement_rules[hic$V1], hic$V1)
hic$V4 <- ifelse(hic$V4 %in% names(replacement_rules), replacement_rules[hic$V4], hic$V4)

#Read chromosome lengths
chromosome_lengths <- read_delim("/path/to/chromosome_lengths_file.genome", 
                                 delim = "\t", escape_double = FALSE, 
                                 col_names = FALSE, trim_ws = TRUE)

#Convert to named vector
chromosome_lengths <- setNames(chromosome_lengths$X2, chromosome_lengths$X1)

#Apply the replacement rules to chromosome lengths
names(chromosome_lengths) <- ifelse(names(chromosome_lengths) %in% names(replacement_rules), replacement_rules[names(chromosome_lengths)], names(chromosome_lengths))

#Calculate the number of bins for each chromosome
chromosome_bins <- sapply(chromosome_lengths, function(length) ceiling(length / bin_size))

#Filter chromosome_bins to include only chromosomes 1-26
chromosome_bins <- chromosome_bins[names(chromosome_bins) %in% as.character(1:26)]

#Generate all possible bin pairs
all_pairs <- expand.grid(chr1 = names(chromosome_bins), bin1 = 1:max(chromosome_bins), 
                         chr2 = names(chromosome_bins), bin2 = 1:max(chromosome_bins))

#Filter out invalid pairs (bins that exceed chromosome length)
all_pairs <- all_pairs %>%
  filter(bin1 <= chromosome_bins[chr1] & bin2 <= chromosome_bins[chr2])

#Convert bin positions to genomic coordinates (0-based)
all_pairs <- all_pairs %>%
  mutate(start1 = (bin1 - 1) * bin_size, end1 = bin1 * bin_size,
         start2 = (bin2 - 1) * bin_size, end2 = bin2 * bin_size)

#Create a data frame with zero interactions for all pairs
all_pairs$int <- 0

#Merge with the existing Hi-C data
hic <- hic %>%
  rename(chr1 = V1, start1 = V2, end1 = V3, chr2 = V4, start2 = V5, end2 = V6, int = V7)

merged_data <- merge(all_pairs, hic, by = c("chr1", "start1", "end1", "chr2", "start2", "end2"), all.x = TRUE)

#Filter merged_data to remove the columns bin1, bin2 and int.x
merged_data <- merged_data %>%
  select(-bin1, -bin2, -int.x)

#Change name of int.y to int
names(merged_data)[names(merged_data) == "int.y"] <- "int"

#Replace NA values with 0
merged_data$int[is.na(merged_data$int)] <- 0

#Convert to GInteractions
convertToGI <- function(df) {
  row.regions <- GRanges(df$chr1, IRanges(df$start1, df$end1)) #interaction start
  col.regions <- GRanges(df$chr2, IRanges(df$start2, df$end2)) #interaction end
  gi <- GInteractions(row.regions, col.regions)
  gi$norm.freq <- df$int #Interaction frequencies
  return(gi)
}

hic.gi <- convertToGI(merged_data)

#Filter data
intra <- filter(merged_data, chr1 == chr2)
inter <- filter(merged_data, chr1 != chr2)

#Add specie column
inter$specie <- "Platypus Inter"
intra$specie <- "Platypus Intra"

#Select interactions for each chromosome
select_interactions <- function(data, col) {
  do.call(rbind, lapply(1:31, function(i) subset(data, data[[col]] == i)))
}

inter <- select_interactions(inter, "chr2")
intra <- select_interactions(intra, "chr1")

#Filter out chromosomes not in 1:26
intra <- intra[intra$chr2 %in% 1:26,]
inter <- inter[inter$chr2 %in% 1:26,]

#Create interaction matrix
matrix <- matrix(, nrow = 26, ncol = 26)

for (i in 1:26) {
  chr_intra <- subset(intra, chr1 == i)
  chr_intra_s <- chr_intra[, c("chr1", "chr2", "int")]
  colnames(chr_intra_s) <- c("chr", "chr2", "int")
  
  chr_inter_1 <- subset(inter, chr1 == i)
  chr_inter_1_s <- chr_inter_1[, c("chr1", "chr2", "int")]
  colnames(chr_inter_1_s) <- c("chr", "chr2", "int")
  
  chr_inter_2 <- subset(inter, chr2 == i)
  chr_inter_2_s <- chr_inter_2[, c("chr2", "chr1", "int")]
  colnames(chr_inter_2_s) <- c("chr", "chr2", "int")
  
  total <- rbind(chr_intra_s, chr_inter_1_s, chr_inter_2_s)
  final <- aggregate(int ~ chr + chr2, data = total, FUN = mean)
  ints <- final[!duplicated(final$chr2),]
  ints$chr2 <- as.numeric(ints$chr2)
  line <- ints[order(ints$chr2),]
  
  #Check if lengths match
  if (length(line$int) != ncol(matrix)) {
    print(paste("Mismatch at chromosome:", i))
    print(paste("Length of line$int:", length(line$int)))
    print(paste("Number of columns in matrix:", ncol(matrix)))
    next
  }
  
  matrix[i,] <- line$int
}

matrix1 <- matrix[1:26, 1:26]

#Get upper triangle
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

upper_tri <- get_upper_tri(matrix1)

#Transform matrix into long format
longData <- melt(upper_tri, na.rm = TRUE)
longData$value <- as.numeric(longData$value)

#Filter longData to remove rows where Var1 is equal to Var2
filtered_longData <- filter(longData, Var1 != Var2)

#Calculate lim1 and lim2
min_value <- min(filtered_longData$value)
max_value <- max(filtered_longData$value)

#Auto limits
lim1 <- floor(min_value / 0.0005) * 0.0005
lim2 <- ceiling(max_value / 0.0005) * 0.0005

#Overriding with manual limits
lim1 <- 0.0185
lim2 <- 0.153

midp <- lim1 + ((lim2 - lim1) / 2)

#Plot data
plot <- ggplot(longData, aes(Var2, Var1, fill = value)) + 
  theme_few() +
  ggtitle(paste("Inter-chromosomal interactions for male platypus (", binsize, " bins)", sep = "")) +
  geom_tile(colour = "grey", size = 0.5) +
  scale_fill_gradient2(name = expression("Mean\nInteractions"), low = "#6570BF", high = "#FF4500", mid = "white", midpoint = midp, limits = c(lim1, lim2)) +
  scale_x_discrete(name = "Chromosomes", limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X1", "X2", "X3", "X4", "X5")) + 
  scale_y_discrete(name = "Chromosomes", limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X1", "X2", "X3", "X4", "X5"))

plot

#Save plot
ggsave(paste0("plots/Platypus_Male_", binsize, "_medblueorange_lims_", lim1, "_", lim2, "_zerowindows.pdf"), plot = plot, width = 10, height = 8, units = "in", dpi = 800)

#Save longData
write.table(longData, file = paste0("data/longDatamale_", binsize, "_zerowindows.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

#Save hic
fwrite(merged_data, file = paste0("data/hicmale_", binsize, "_zerowindows.txt"), quote = FALSE, sep = "\t")