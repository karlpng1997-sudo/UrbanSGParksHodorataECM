
# Set path

setwd("C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025")

# Load the CSV file
df <- read.csv("OTU_table_GenusSpecies.csv", check.names = FALSE)

# Function to rename OTU columns
rename_otu <- function(old_name) {
  # Extract OTU number
  otu_num <- sub(".*OTU\\s*_?\\s*(\\d+).*", "\\1", old_name)
  
  # Extract genus and species
  genus <- if (grepl("g__", old_name)) sub(".*g__([A-Za-z0-9]+).*", "\\1", old_name) else NA
  species <- if (grepl("s__", old_name)) sub(".*s__([A-Za-z0-9]+).*", "\\1", old_name) else NA
  
  # Construct new name
  if (!is.na(genus) & !is.na(species)) {
    new_name <- paste0(genus, "_", species, "_OTU", otu_num)
  } else if (!is.na(genus)) {
    new_name <- paste0(genus, "_OTU", otu_num)
  } else {
    new_name <- paste0("NA_OTU", otu_num)
  }
  
  return(new_name)
}

# Apply the function to all column names
colnames(df) <- sapply(colnames(df), rename_otu)

# Save the renamed table
write.csv(df, "OTU_table_GenusSpecies_renamed_old.csv", row.names = FALSE, quote = FALSE)

# 1) Load the table (don't convert names)
df <- read.csv("OTU_table_GenusSpecies_renamed.csv", check.names = FALSE, stringsAsFactors = FALSE)

orig_names <- colnames(df)

# Extract genus as the first chunk before the first underscore
genus_vec <- sub("^([^_]+).*", "\\1", orig_names)

# named numeric vector for counters (starts empty)
counts <- numeric(0)

new_names <- character(length(orig_names))

for (i in seq_along(orig_names)) {
  g <- genus_vec[i]
  # initialize or increment
  if (is.na(counts[g])) {
    counts[g] <- 1
  } else {
    counts[g] <- counts[g] + 1
  }
  # robustly replace the OTU number at end (handles "OTU3", "OTU_3", "OTU _ 3", "OTU 3")
  if (grepl("OTU", orig_names[i], ignore.case = FALSE)) {
    new_names[i] <- sub("OTU\\s*[_ ]*\\s*\\d+\\s*$", paste0("OTU", counts[g]), orig_names[i])
  } else {
    # if no OTU token was present, append new OTU number
    new_names[i] <- paste0(orig_names[i], "_OTU", counts[g])
  }
}

# Assign and save
colnames(df) <- new_names
write.csv(df, "OTU_table_GenusSpecies_renumbered_by_genus_order_old.csv", row.names = FALSE, quote = FALSE)

df <- read.csv("OTU_table_GenusSpecies_renumbered_by_genus_order.csv", row.names = 1)

df_prop <- df / rowSums(df)

write.csv(df_prop, "OTU_table_GenusSpecies_proportions.csv", quote = FALSE)

# Load your proportion or count table
df <- read.csv("OTU_table_GenusSpecies_proportions.csv", check.names = FALSE, row.names = 1)

# Define ECM genera
ecm_genera <- c("Russula", "Tomentella", "Lactarius", "Inocybe", "Lactifluus", "Clavulina", "Scleroderma")

# Identify ECM vs non-ECM columns
is_ecm <- sapply(colnames(df), function(x) any(sapply(ecm_genera, function(g) grepl(g, x, ignore.case = TRUE))))

# Create a new column summing all non-ECM columns
df$`Non-ECM fungi` <- rowSums(df[, !is_ecm, drop = FALSE])

# Keep only ECM columns + the new one
df_ecm_summary <- cbind(df[, is_ecm, drop = FALSE], df[, "Non-ECM fungi", drop = FALSE])

# Save the result
write.csv(df_ecm_summary, "OTU_table_ECM_and_nonECM_old.csv", quote = FALSE)

# Open libraries

library(ggplot2)
library(tidyr)
library(dplyr)

# Read data
df <- read.csv("OTU_table_ECM_and_nonECM.csv", row.names = 1, check.names = FALSE)

# Convert row names (samples) to a column
df$Sample <- rownames(df)

# Convert from wide to long format
df_long <- pivot_longer(df, 
                        cols = -Sample, 
                        names_to = "OTU", 
                        values_to = "Proportion")

group_colors <- c('Clavulina_OTU01' = 'red',
                  'Inocybe_OTU01' = 'purple',
                  'Lactarius_OTU01' = 'palegreen1',
                  'Lactifluus_OTU01' = 'blue',
                  'Russula_OTU01' = 'pink1',
                  'Russula_OTU02' = 'palevioletred2',
                  'Russula_OTU03' = 'pink3',
                  'Russula_OTU04_violeipes' = 'palevioletred4',
                  'Scleroderma_OTU01' = 'gold1',
                  'Scleroderma_OTU02' = 'gold3',
                  'Tomentella_OTU01' = 'grey3',
                  'Tomentella_OTU02' = 'grey6',
                  'Tomentella_OTU03' = 'grey7',
                  'Tomentella_OTU04' = 'grey14',
                  'Tomentella_OTU05' = 'grey16',
                  'Tomentella_OTU06' = 'grey18',
                  'Tomentella_OTU07' = 'grey20',
                  'Tomentella_OTU08' = 'grey23',
                  'Tomentella_OTU09' = 'grey25',
                  'Tomentella_OTU10' = 'grey27',
                  'Tomentella_OTU11' = 'grey29',
                  'Tomentella_OTU12' = 'grey31',
                  'Tomentella_OTU13' = 'grey33',
                  'Tomentella_OTU14_guineensis' = 'grey70',
                  'Tomentella_OTU15_parmastoana' = 'grey80',
                  '*Non-ECM fungi' = 'beige')

ggplot(df_long, aes(x = Sample, y = Proportion, fill = OTU)) +
  geom_bar(stat = "identity") +
  theme(
    legend.position = "right"
  ) +
  labs(
    x = "Park Location Sample",
    y = "Relative Abundance"
  )+
  scale_fill_manual(values = group_colors)+
  guides(fill = guide_legend(title = "Fungal OTUs", ncol = 3, order = 1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


df_long <- df_long %>%
  mutate(Genus = sub("_.*", "", OTU))  # Extract text before first underscore

ggplot(df_long, aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "Proportion of Fungal Genera per Sample",
    x = "Sample",
    y = "Proportion"
  )
