#Set working directory
setwd("/path/to/working/dir")

library(dplyr)
library(RColorBrewer)
library(ggplot2)

#Define binsize variable
binsize <- "500kb"

#Read in male data
Male_subset <- read.table(paste0("data/longDatamale_", binsize, "_zerowindows.txt"), header = TRUE)

#Read in female data
Female <- read.table(paste0("data/longDatafemale_", binsize, "_zerowindows.txt"), header = TRUE)

#Remove all rows in Male_subset where Var1 is equal to Var2
Male_subset <- Male_subset[Male_subset$Var1 != Male_subset$Var2, ]

#Remove all rows in Female where Var1 is equal to Var2
Female <- Female[Female$Var1 != Female$Var2, ]

#Sort the dataframe by the 'value' column in descending order and select the top 10 rows
top_10_Male_subset <- Male_subset %>%
  arrange(desc(value)) %>%
  slice_head(n = 10)

#Sort the dataframe by the 'value' column in descending order and select the top 10 rows
top_10_Female <- Female %>%
  arrange(desc(value)) %>%
  slice_head(n = 10)

#Merge the two data frames on Var1 and Var2
common_rows <- merge(top_10_Male_subset, top_10_Female, by = c("Var1", "Var2"))

#Check if there are any common rows
if (nrow(common_rows) > 0) {
  print("There are common rows between the two data frames:")
  print(common_rows)
} else {
  print("There are no common rows between the two data frames.")
}

#Recode Var1 and Var2
top_10_Female <- top_10_Female %>%
  mutate(
    Var1 = as.character(Var1),
    Var2 = as.character(Var2),
    Var1 = case_when(
      Var1 == "22" ~ "X1",
      Var1 == "23" ~ "X2",
      Var1 == "24" ~ "X3",
      Var1 == "25" ~ "X4",
      Var1 == "26" ~ "X5",
      TRUE ~ Var1  #Keep original if no match
    ),
    Var2 = case_when(
      Var2 == "22" ~ "X1",
      Var2 == "23" ~ "X2",
      Var2 == "24" ~ "X3",
      Var2 == "25" ~ "X4",
      Var2 == "26" ~ "X5",
      TRUE ~ Var2  #Keep original if no match
    )
  )

top_10_Male_subset <- top_10_Male_subset %>%
  mutate(
    Var1 = as.character(Var1),
    Var2 = as.character(Var2),
    Var1 = case_when(
      Var1 == "22" ~ "X1",
      Var1 == "23" ~ "X2",
      Var1 == "24" ~ "X3",
      Var1 == "25" ~ "X4",
      Var1 == "26" ~ "X5",
      TRUE ~ Var1  #Keep original if no match
    ),
    Var2 = case_when(
      Var2 == "22" ~ "X1",
      Var2 == "23" ~ "X2",
      Var2 == "24" ~ "X3",
      Var2 == "25" ~ "X4",
      Var2 == "26" ~ "X5",
      TRUE ~ Var2  #Keep original if no match
    )
  )

#Combine Var1 and Var2 to form unique words
top_10_Female$Word <- paste(top_10_Female$Var1, "and", top_10_Female$Var2, sep = " ")

#Combine Var1 and Var2 for Male
top_10_Male_subset$Word <- paste(top_10_Male_subset$Var1, "and", top_10_Male_subset$Var2, sep = " ")

#Assign ranks (based on 'value') and identify the list (Female/Male)
top_10_Female <- top_10_Female %>% mutate(Rank = rank(desc(value)), List = "Female")
top_10_Male_subset <- top_10_Male_subset %>% mutate(Rank = rank(desc(value)), List = "Male subset")

#Combine the data frames
combined_df <- bind_rows(top_10_Female, top_10_Male_subset)

#Assign unique colors to each word using custom HEX codes
unique_words <- unique(combined_df$Word)
custom_hex_colors <- c("#CAE1FF", "#FF4500", "#6570BF", "#FFDAB9", "#000080", "#FF8F5C")
color_palette <- colorRampPalette(custom_hex_colors)(length(unique_words))
colors <- setNames(color_palette, unique_words)

#Plotting
ggplot() +
  geom_point(data = top_10_Female, aes(x = List, y = Rank, color = Word), size = 5) +
  geom_point(data = top_10_Male_subset, aes(x = List, y = Rank, color = Word), size = 5) +
  geom_line(data = combined_df %>% filter(Word %in% intersect(top_10_Female$Word, top_10_Male_subset$Word)),
            aes(x = List, y = Rank, group = Word, color = Word), size = 1) +
  geom_text(data = combined_df, aes(x = List, y = Rank, label = Word, color = Word,
                                    hjust = ifelse(List == "Female", 1.2, -0.2)),
            size = 3.5, vjust = 0.5) +  #Adjust hjust and vjust to position the labels
  scale_color_manual(values = colors) +
  scale_y_reverse() +
  labs(title = "Top 10 non-self interchromosomal interactions",
       y = "Chromosome pair interaction strength",
       x = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),  #Remove y-axis text
    axis.ticks.y = element_blank(),  #Remove y-axis ticks
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank()   #Remove minor gridlines
  )

#Save plot
ggsave(paste0("plots/Top_10_non-self_interchromosomal_interactions_normalised_", binsize, "_zerowindows_blueorange.pdf"), width = 4.5, height = 3.5, dpi = 600)
