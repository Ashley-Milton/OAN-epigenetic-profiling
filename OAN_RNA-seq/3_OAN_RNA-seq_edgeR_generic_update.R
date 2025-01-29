#Set working directory
setwd("/path/to/working/dir")

#Load packages
library(edgeR)
library(tidyverse)
library(readr)
library(stringr)
library(gplots)

#Read in Subread featureCounts data
counts <- read.table("featureCounts_file.txt", header=T)

#Select relevant columns
counts <- counts %>% select(-Start, -End)

#Read in genome annotation file
GeneGFF <- read_delim("/path/to/genome/directory/annotation_file.gtf", 
                                                                                     delim = "\t", escape_double = FALSE, 
                                                                                     col_names = FALSE, trim_ws = TRUE, skip = 4)

GeneGFF <- GeneGFF %>% filter(!grepl("##",X1)) #For removing unneeded last line

#Filter for genes only and select relevant columns
GeneNames <- GeneGFF %>% filter(X3 =="gene") %>% select(X1,X4,X9)

#Update column names
colnames(GeneNames) <- c("Chromosome","Start","ID") 

#Separate the ID column into categories using key-value pairs
GeneNames_separated <- GeneNames %>%
  separate_rows(ID, sep = "; ") %>%
  separate(ID, into = c("key", "value"), sep = " ", extra = "merge") %>%
  mutate(key = trimws(key), value = trimws(gsub('"', '', value))) %>%
  pivot_wider(names_from = key, values_from = value, values_fn = list(value = list)) %>%
  mutate(across(where(is.list), ~ sapply(., function(x) if (is.null(x)) NA else paste(unlist(x), collapse = "; "))))

#Select relevant columns
GeneRef <- GeneNames_separated %>% select(gene_id, gene)

#Perform a left join of counts and GeneRef to add gene column from GeneRef
oanCounts <- left_join(counts, GeneRef, by = c("Geneid" = "gene_id"))

#Set row names to the gene column
rownames(oanCounts) <- make.names(oanCounts$gene, unique = TRUE)

#Select the relevant columns (excluding the gene column)
oanCounts <- oanCounts %>% select(-gene)

#Rename specific columns by position
colnames(oanCounts)[c(5, 6)] <- c("Female", "Male")

#Select only the numeric columns (Female and Male)
count_matrix <- oanCounts %>% select(Female, Male)

#Create a DGEList object
d <- DGEList(counts = as.matrix(count_matrix), group = c("Female", "Male"))

#View the DGEList object
print(head(d))

#Normalize the counts
d <- calcNormFactors(d)

#Create a dataframe from the normalized counts
df1 <- data.frame(d$counts)

#View quantiles of row sums
print(quantile(rowSums(df1)))

#Convert the dataframe to a matrix
data_matrix <- data.matrix(df1)

#Calculate CPM (Counts Per Million)
cpm_values <- cpm(d)

#Convert CPM values to a dataframe again for easier manipulation
cpm_df <- as.data.frame(cpm_values)

#Define gene(s) of interest
genes_of_interest <- "AMH"

#Extract CPM values for the genes of interest
expression_comparison <- cpm_df[rownames(cpm_df) %in% genes_of_interest, ]

#View the M/F expression levels (CPM) of the gene of interest
print(expression_comparison)

#Convert the data to long format for ggplot2
expression_long <- expression_comparison %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")

#Visualize the expression levels using a bar plot
expression_long %>%
  filter(gene == "AMH") %>%
  ggplot(aes(x = gene, y = expression, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(sample, "\n", sprintf("%.2f", expression))), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3) +
  labs(x = NULL, y = "Counts Per Million (CPM)") +
  theme_classic() +
  ylim(c(0, 6)) +
  scale_fill_manual(values = c("Female" = "#000080", "Male" = "#FF4500"), name = NULL) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 14))

ggsave("/path/to/outdir/CPM_plot.pdf", dpi=800, width = 5, height = 5)
