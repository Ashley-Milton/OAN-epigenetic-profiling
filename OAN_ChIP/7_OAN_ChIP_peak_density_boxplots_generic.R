#Set the working directory
setwd("/path/to/peak/files")

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(conover.test)
library(stringr)

#Define headers for broadPeak and narrowPeak files
broadPeak_headers <- c("chromosome", "start_coordinate", "end_coordinate", "name", "score", "strand", "signalValue", "pValue", "qValue")
narrowPeak_headers <- c("chromosome", "start_coordinate", "end_coordinate", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

#The following lines deal with specific file naming patterns

#List all .broadPeak and .narrowPeak files
all_files <- list.files(pattern = "\\.(broadPeak|narrowPeak)$")

#Define the mapping of input file names to shorter, standard names
name_mapping <- list(
  "F_H3K27me3_peaks.broadPeak" = "F_H3K27me3_ours",
  "F_H3K4me_german_peaks.broadPeak" = "F_H3K4me1",
  "F_H3K9me2_peaks.broadPeak" = "F_H3K9me2",
  "F_H3K9me3_peaks.broadPeak" = "F_H3K9me3",
  "F_H4K20me1_peaks.broadPeak" = "F_H4K20me1",
  "K27me3_F_366_broad_peaks.broadPeak" = "F_H3K27me3",
  "K27me3_M_366_broad_peaks.broadPeak" = "M_H3K27me3",
  "M_H3K27me3_peaks.broadPeak" = "M_H3K27me3_ours",
  "M_H3K4me_german_peaks.broadPeak" = "M_H3K4me1",
  "M_H3K9me2_peaks.broadPeak" = "M_H3K9me2",
  "M_H3K9me3_peaks.broadPeak" = "M_H3K9me3",
  "M_H4K20me1_peaks.broadPeak" = "M_H4K20me1",
  "F_H3K27ac_german_peaks.narrowPeak" = "F_H3K27ac",
  "F_H3K4me3_german_peaks.narrowPeak" = "F_H3K4me3",
  "M_H3K27ac_german_peaks.narrowPeak" = "M_H3K27ac",
  "M_H3K4me3_german_peaks.narrowPeak" = "M_H3K4me3"
)

#Read each file into a separate data frame
for (file in all_files) {
  #Get the new name from the mapping
  var_name <- name_mapping[[file]]
  
  #Read the file with appropriate headers
  if (grepl("\\.broadPeak$", file)) {
    df <- fread(file, col.names = broadPeak_headers)
  } else {
    df <- fread(file, col.names = narrowPeak_headers)
  }
  
  #Filter to only keep qValues above 1.3
  df <- df %>% filter(qValue > 1.3)
  
  #Assign the data frame to a variable with the new name
  assign(var_name, df)
}

#Remove unwanted dataframes
rm(df, M_H3K27me3_ours, F_H3K27me3_ours)

#Load the chromosome name lookup table
lookup_table <- fread("/path/to/lookup_table.csv", col.names = c("old_name", "new_name"))

#Load the chromosome lengths file
chrom_lengths <- fread("/path/to/chromsizes.genome", col.names = c("chromosome", "length"))

#Function to map chromosome names using the lookup table
map_chromosome_names <- function(chromosome) {
  mapped_name <- lookup_table[lookup_table$old_name == chromosome, "new_name"]
  if (nrow(mapped_name) == 0) {
    return(chromosome)
  } else {
    return(mapped_name$new_name)
  }
}

#Function to create 1 Mb windows for a given chromosome with 0-based indexing
create_windows <- function(chromosome, length) {
  GRanges(seqnames = chromosome,
          ranges = IRanges(start = seq(0, length - 1, by = 1e6),
                           end = pmin(seq(1e6, length + 1e6, by = 1e6), length)))
}

#Create a common set of 1 Mb windows for each chromosome based on the chromosome lengths file
windows <- lapply(1:nrow(chrom_lengths), function(i) {
  chr <- chrom_lengths$chromosome[i]
  length <- chrom_lengths$length[i]
  create_windows(chr, length)
})
all_windows <- do.call(c, windows)

#Function to determine if a row is X-specific (not in a PAR)
is_Xspec <- function(chr, start, end) {
  (chr == "X1" & start > 45100000) |  #X1
    (chr == "X2" & start > 5700000 & end < 20000000) |  #X2
    (chr == "X3" & start > 18500000 & end < 31900000) |  #X3
    (chr == "X4" & start > 5300000) |  #X4
    (chr == "X5" & start > 8500000 & end < 69200000)  #X5
}

#Function to process each histone modification
process_histone_mod <- function(histone_mod) {
  histone_mod <- histone_mod %>%
    mutate(chromosome = sapply(chromosome, map_chromosome_names)) %>%
    mutate(chr_type = case_when(
      chromosome %in% paste0("X", 1:5) ~ "Xs",
      chromosome %in% paste0("Y", 1:5) ~ "Ys",
      chromosome == "Mitochondrial" ~ "mt",
      chromosome %in% as.character(1:21) ~ "autosomes",
      TRUE ~ "unknown"
    )) %>%
    filter(!(chr_type %in% c("Ys", "unknown", "mt")))
  
  peaks <- GRanges(seqnames = histone_mod$chromosome,
                   ranges = IRanges(start = histone_mod$start_coordinate,
                                    end = histone_mod$end_coordinate))
  
  counts <- countOverlaps(all_windows, peaks)
  
  autosomes_windows <- all_windows[seqnames(all_windows) %in% as.character(1:21)]
  autosomes_counts <- countOverlaps(autosomes_windows, peaks)
  autosomes_windows_count <- length(autosomes_windows)
  
  Xs_windows <- all_windows[seqnames(all_windows) %in% paste0("X", 1:5)]
  
  Xspec_windows <- Xs_windows[is_Xspec(seqnames(Xs_windows), start(Xs_windows), end(Xs_windows))]
  Xspec_counts <- countOverlaps(Xspec_windows, peaks)
  Xspec_windows_count <- length(Xspec_windows)
  
  PARs_windows <- Xs_windows[!is_Xspec(seqnames(Xs_windows), start(Xs_windows), end(Xs_windows))]
  PARs_counts <- countOverlaps(PARs_windows, peaks)
  PARs_windows_count <- length(PARs_windows)
  
  list(
    overall_windows = length(all_windows),
    autosomes_windows_count = autosomes_windows_count,
    Xspec_windows_count = Xspec_windows_count,
    PARs_windows_count = PARs_windows_count,
    autosomes_counts = autosomes_counts,
    Xspec_counts = Xspec_counts,
    PARs_counts = PARs_counts
  )
}

#Initialise lists to store counts
autosomes_counts_list <- list()
Xspec_counts_list <- list()
PARs_counts_list <- list()

#Process each histone modification
results <- list()
for (name in names(name_mapping)) {
  var_name <- name_mapping[[name]]
  if (exists(var_name)) {
    histone_mod <- get(var_name)
    result <- process_histone_mod(histone_mod)
    results[[var_name]] <- result
    
    #Save autosomes_counts, Xspec_counts, and PARs_counts as data frames
    autosomes_counts_list[[var_name]] <- data.frame(window = seq_along(result$autosomes_counts), counts = result$autosomes_counts)
    Xspec_counts_list[[var_name]] <- data.frame(window = seq_along(result$Xspec_counts), counts = result$Xspec_counts)
    PARs_counts_list[[var_name]] <- data.frame(window = seq_along(result$PARs_counts), counts = result$PARs_counts)
  }
}

#Combine autosomes counts into a single data frame
autosomes_df <- do.call(rbind, lapply(names(autosomes_counts_list), function(name) {
  df <- autosomes_counts_list[[name]]
  df$histone_mod <- name
  df$type <- "autosomes"
  df
}))

#Combine Xs counts into a single data frame
Xspec_df <- do.call(rbind, lapply(names(Xspec_counts_list), function(name) {
  df <- Xspec_counts_list[[name]]
  df$histone_mod <- name
  df$type <- "Xs"
  df
}))

#Combine PARs counts into a single data frame
PARs_df <- do.call(rbind, lapply(names(PARs_counts_list), function(name) {
  df <- PARs_counts_list[[name]]
  df$histone_mod <- name
  df$type <- "PARs"
  df
}))

#Combine all data frames into one
combined_df <- rbind(autosomes_df, Xspec_df, PARs_df)

#Separate the histone_mod column into sex and histone_mod columns
combined_df <- combined_df %>%
  separate(histone_mod, into = c("sex", "histone_mod"), sep = "_", extra = "merge")

#Create a new column to combine sex and type to use in plotting
combined_df <- combined_df %>%
  mutate(group = paste(sex, type, sep = "_"))

#Create colour palette
colors <- c(colorRampPalette(c("lightsteelblue1", "navy"))(3), colorRampPalette(c("peachpuff", "orangered1"))(3))

#Define the desired order for the histone_mod facets
desired_order <- c("H3K4me1","H3K4me3", "H3K27ac", "H4K20me1", "H3K9me2", "H3K9me3", "H3K27me3")

#Convert histone_mod to a factor with the desired order
combined_df$histone_mod <- factor(combined_df$histone_mod, levels = desired_order)

#Calculate medians for each group
medians <- combined_df %>%
  group_by(sex, histone_mod, type, group) %>%
  summarize(median_count = median(counts), .groups = 'drop')

#Calculate means for each group
means <- combined_df %>%
  group_by(sex, histone_mod, type, group) %>%
  summarize(mean_count = mean(counts), .groups = 'drop')

#Create the boxplot with free y-scales and faceting by histone_mod and sex
ggplot(combined_df, aes(x = type, y = counts, fill = group)) +
  geom_boxplot(notch = T, outlier.size = 0.1) +
  stat_summary(fun = mean, geom = "point", color = "black") +
  geom_text(data = means, aes(x = type, y = mean_count, label = signif(mean_count, 2)), 
            position = position_dodge(width = 1), vjust = -6.5, hjust = -0.35, size = 3) +
  labs(title = "Boxplot of Peaks/Mb for Each Histone Modification",
       x = NULL,
       y = "Peaks/Mb",
       fill = NULL) +
  scale_fill_manual(values = colors, 
                    labels = c("F_autosomes" = "Female autosomes", 
                               "F_Xs" = "Female X-specific", 
                               "F_PARs" = "Female PARs",
                               "M_autosomes" = "Male autosomes", 
                               "M_Xs" = "Male X-specific",
                               "M_PARs" = "Male PARs")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 8),
        panel.grid = element_blank(),  #Remove gridlines
        strip.background = element_rect(fill = NA, color = "black")) +
  facet_grid(histone_mod ~ sex, scales = "free_y", labeller = labeller(sex = c(F = "Female", M = "Male")))

ggsave("/path/to/outdir/Peaks_per_Mb_boxplots.pdf", dpi = 300, width = 16, height = 24, units = "cm")


fwrite(means, file = "/path/to/outdir/means.txt", sep = "\t")


##Stats tests

combined_df <- combined_df %>%
  mutate(histone_group = paste(histone_mod, group, sep = "_"))

fwrite(combined_df, "/path/to/outdir/combined_df.txt", sep = "\t")

#Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(counts ~ histone_group, data = combined_df)
print(kruskal_test_result)

#If the Kruskal-Wallis test is significant, proceed with Conover-Iman test
if (kruskal_test_result$p.value < 0.05) {
  conover_test_result <- conover.test(combined_df$counts, combined_df$histone_group, method = "bonferroni")
  
  #Extract the adjusted p-values and comparisons
  adjusted_p_values <- conover_test_result$P.adjusted
  comparisons <- conover_test_result$comparisons
  
  #Combine into a data frame
  results_df <- data.frame(comparisons, adjusted_p_values)
  
  #Filter for significant comparisons (adjusted p-value < 0.05)
  significant_results <- subset(results_df, adjusted_p_values < 0.05)
  
  #Filter for non significant comparisons (adjusted p-value < 0.05)
  non_significant_results <- subset(results_df, adjusted_p_values >= 0.05)
  
  #Print significant results
  print(significant_results)
} else {
  cat("No significant differences found by Kruskal-Wallis test.\n")
}

fwrite(results_df, "/path/to/outdir/Conover_results_df.txt", sep = "\t")

#Filter the data for within-histone comparisons only
filtered_df <- results_df %>%
  filter(
    str_extract(comparisons, "^[^_]+_[^_]+") == str_extract(comparisons, "(?<= - )[^_]+_[^_]+")
  )
