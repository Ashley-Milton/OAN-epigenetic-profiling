#Set the working directory
setwd("/path/to/indir")

library(tidyverse)
library(data.table)
library(scales)
library(readr)
library(grid)

#Define histone variable
histones <- c("H3K4me3_german","H3K4me_german","H3K9me2","H3K9me3","H3K27ac_german","H3K27me3","H4K20me1")

#Load the .gtf file and extract the positions for the AMHY gene
gtf <- read_delim("/path/to/genome/refgenome.gtf", 
                  delim = "\t", escape_double = FALSE, 
                  col_names = FALSE, trim_ws = TRUE, skip = 4)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#Filter for gene features and extract the AMHY gene
amh_gene <- gtf %>% filter(feature == "gene" & grepl("AMH;", attribute))

#Extract the start and end positions of the AMH gene
amh_start <- min(amh_gene$end)
amh_end <- max(amh_gene$start)
amh_chrom <- unique(amh_gene$seqname)

#Define the region around the AMH TSS (10 kb upstream and downstream)
region_start <- amh_start + 10000
region_end <- amh_end - 1000

#Load the lookup table
lookup_table <- fread("/path/to/lookup_table.csv", col.names = c("old_name", "new_name"))

for (histone in histones) {
  data_M <- fread(paste0("/path/to/macs2/output/M_", histone, "_treat_pileup.bdg"), header = FALSE)
  
  #Rename columns to match the expected format
  colnames(data_M) <- c("chromosome", "start", "end", "read_depth")
  
  #Join with lookup table to rename chromosomes
  data_M <- data_M %>% 
    left_join(lookup_table, by = c("chromosome" = "old_name")) %>% 
    mutate(chromosome = ifelse(is.na(new_name), chromosome, new_name)) %>% 
    select(-new_name)
  
  #Filter for the chromosome of interest
  data_M <- data_M %>% filter(chromosome == amh_chrom)
  
  #Convert start and end columns to numeric
  data_M <- data_M %>% mutate(across(c(start, end), as.numeric))
  
  #TSS plots for AMHY
  data_M <- data_M %>% mutate(mod = histone, sex = "Male")
  
  data_M %>% 
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = read_depth), fill = "#FF8F5C") +
    geom_vline(xintercept = amh_start, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = amh_end, linetype = "dashed", color = "black", linewidth = 0.5) +
    annotate("text", x = amh_start, y = Inf, label = "AMHY TSS", vjust = 2, hjust = 1, angle = 90, color = "black", size = 4) +
    annotate("text", x = amh_end, y = Inf, label = "AMHY end", vjust = 2, hjust = 1, angle = 90, color = "black", size = 4) +
    annotate("segment", x = amh_start, xend = amh_end, y = 0, yend = 0, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "black") +
    labs(x = "Position on Y5 around TSS of AMHY (kb)", y = paste("Read depth of", histone, "(male)")) +
    coord_cartesian(x = c(region_end, region_start), y = c(0, 20))
    theme_classic() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          aspect.ratio = 1,
          legend.title = element_blank())
  
  ggsave(paste("path/to/outdir/TSS_AMHY_M_", histone, "_plot.pdf", sep = ""), width = 7, height = 7, dpi = 800)
  
}