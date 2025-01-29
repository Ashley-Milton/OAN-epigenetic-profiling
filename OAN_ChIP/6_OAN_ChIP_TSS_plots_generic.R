#Set working directory and load packages
setwd("/path/to/plotting/data")

library(tidyverse)
library(data.table)
library(scales)
library(fuzzyjoin)

#Define histone variable
histones <- c("H3K27ac", "H3K27me3", "H3K4me3", "H3K4me", "H3K9me2", "H3K9me3", "H4K20me1")

#Load common data before loop

#Read the file containing the number of mapped reads
mapped_reads <- read.csv("mapped_read_counts.csv")

#Read the file containing the PARs
PARs <- read.csv("PAR_boundaries.csv")

#Read the file containing the chromosome name lookup table
lookup <- read.csv("Platypus_lookup_table.csv")

#Loop through histones
for (histone in histones) {

  #Read data in
  data_F <- fread(paste0("F_", histone, "_info_simplified_realpos.tsv"), select = c("position", "relative_pos", "height"))
  data_M <- fread(paste0("M_", histone, "_info_simplified_realpos.tsv"), select = c("position", "relative_pos", "height"))

  #Normalising
  
  #Get the number of mapped reads for the current histone and sex
  mapped_reads_F <- filter(mapped_reads, histone_type == histone, sex == "Female")$mapped_reads
  mapped_reads_M <- filter(mapped_reads, histone_type == histone, sex == "Male")$mapped_reads
  
  #Calculate the scaling factor
  scaling_factor <- max(mapped_reads_F, mapped_reads_M) / min(mapped_reads_F, mapped_reads_M)
  
  #Determine which dataset has fewer reads
  if (mapped_reads_F < mapped_reads_M) {
    #If female data has fewer reads, scale up the female data
    data_F$height <- data_F$height * scaling_factor
  } else {
    #If male data has fewer reads, scale up the male data
    data_M$height <- data_M$height * scaling_factor
  }
  
  #Further data processing

  process_data <- function(data, lookup) {
    #Convert dataframe to a data.table
    setDT(data)
    
    #Use tstrsplit to split the 'position' column
    data[, c("chr", "start", "end") := tstrsplit(position, ",", fixed=TRUE)]
    
    #Remove the original 'position' column
    data[, position := NULL]
    
    #Change chrom names - merge using lookup table
    data <- data[lookup, on = "chr"]
    
    #Delete old chr column
    data[, chr := NULL]
    
    #Rename the 'chromosome' column to 'chr'
    setnames(data, old = "chromosome", new = "chr")
    
    #Filter data to remove values in the chr column that are "Mitochondrial" or start with "Y"
    data <- data[!chr %in% "Mitochondrial" & !grepl("^Y", chr)]
    
    return(data)
  }

  #FEMALE
  data_F <- process_data(data_F, lookup)

  #MALE
  data_M <- process_data(data_M, lookup)
  
  #Here assign PAR or X-specific depending on whether the peak is in the PAR or not
  
  #Convert your PAR dataframe to data.table
  setDT(PARs)

  assign_PAR_info <- function(data, PARs) {
    #Convert the 'start' and 'end' columns in data and PARs to integer
    data[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
    PARs[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
    
    #Create the new column PAR_info in data that is autosomes by default
    data[, PAR_info := "autosomes"]
    
    #Update PAR_info to "X-specific" for rows where chr starts with "X"
    data[grepl("^X", chr), PAR_info := "X-specific"]
    
    #Perform a non-equi join to update PAR_info to "PAR" where chr values match and start and end positions are within the start and end positions in PARs
    data[PARs, on = .(chr, start >= start, end <= end), PAR_info := "PAR"]
    
    return(data)
  }

  #FEMALE
  data_F <- assign_PAR_info(data_F, PARs)

  #MALE
  data_M <- assign_PAR_info(data_M, PARs)
  
  #Filtering data for TSS plots ---------------------------------------------------------------
  
  #Filtering for at least two datapoints per relative position, and
  #taking median value at each position relative to TSS

  calc_pos_means <- function(data, sex) {
    meandata <- data[, count := .N, by = .(relative_pos, PAR_info)
                    ][count > 2
                      ][, .(mean_height = mean(height)), by = .(relative_pos, PAR_info)
                        ][, sex := sex]
    return(meandata)
  }

  #FEMALE
  meandata_F <- calc_pos_means(data_F, "Female")

  #MALE
  meandata_M <- calc_pos_means(data_M, "Male")
  
  #Putting female and male data together in one dataframe
  meandata <- rbindlist(list(meandata_F, meandata_M))
  meandata[, sex_PAR_info := paste(sex, PAR_info, sep = " ")]

  rm(data_F, data_M, meandata_F, meandata_M)
  gc()
  
  #TSS plot
  
  #Create color palette
  color_palette <- c(colorRampPalette(c("lightsteelblue1", "navy"))(3), colorRampPalette(c("peachpuff", "orangered1"))(3))
  
  #Plot
  meandata %>% ggplot(aes(x=relative_pos, y=mean_height, color=sex_PAR_info)) +
    geom_smooth(aes(fill=sex_PAR_info), se=T, method="loess", span = 0.35) +
    labs(x="Position relative to TSS (kb)", y=paste("Mean peak height (normalised) of", histone)) + 
    coord_cartesian(x=c(-10000,10000))+
    theme_classic()+
    scale_color_manual(values = color_palette)+
    scale_fill_manual(values = color_palette)+
    theme(legend.position="bottom")+
    scale_x_continuous(labels = unit_format(unit="", scale = 1e-3))+
    theme(legend.text = element_text(size = 14),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          aspect.ratio = 1,
          legend.position = "right",
          legend.title = element_blank())
  
  ggsave(paste("/path/to/outdir/TSS_", histone, "_plot.pdf", sep = ""), width = 9, height = 7, dpi = 800)

  rm(meandata)
  gc()

}