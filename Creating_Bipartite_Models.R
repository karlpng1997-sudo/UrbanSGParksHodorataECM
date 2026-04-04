
# Set path

setwd("C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025")

# Load libraries

library(RColorBrewer)
library(bipartite)
library(dplyr)
library(ggplot2)

# Load and prepare data
otu <- read.csv("OTU_table_ECM_and_nonECM_converted_bipartite.csv",
                row.names = 1, check.names = FALSE)

otu <- as.matrix(otu)
mode(otu) <- "numeric"

# You can reorder to group similar sites or OTUs
# For example, by row/column totals (more frequent = central)

otu_sorted <- otu[order(rowSums(otu), decreasing = TRUE),
                  order(colSums(otu), decreasing = TRUE)]

# Basic bipartite plot
plotweb(
  otu_sorted,
  method = "normal",          # standard bipartite layout
  bor.col.interaction = NA,   # remove border on links
  col.low = "skyblue",   # lower level (samples)
  col.high = "tomato",
  text.rot = 90,              # rotate OTU labels for readability
  labsize = 5               # adjust label size if crowded
)

# Load libraries

library(dplyr)
library(readr)

# Read your data
otu <- read_csv("OTU_table_ECM_and_nonECM_converted_locationbipartite.csv")

# Check the structure
head(otu)

# Combine OTU presence/absence by site
otu_by_site <- otu %>%
  group_by(site) %>%
  summarise(across(contains("OTU"), ~ as.numeric(any(. == 1)))) %>%
  ungroup()

# View the result
head(otu_by_site)

write_csv(otu_by_site, "OTU_by_site_presence_absence.csv")

# Load and prepare data
otu <- read.csv("OTU_by_site_presence_absence.csv",
                row.names = 1, check.names = FALSE)

otu <- as.matrix(otu)
mode(otu) <- "numeric"

# You can reorder to group similar sites or OTUs
# For example, by row/column totals (more frequent = central)
otu_sorted <- otu[order(rowSums(otu), decreasing = TRUE),
                  order(colSums(otu), decreasing = TRUE)]

# Basic bipartite plot
plotweb(
  otu_sorted,
  method = "normal",          # standard bipartite layout
  bor.col.interaction = NA,   # remove border on links
  col.low = "skyblue",   # lower level (samples)
  col.high = "tomato",
  text.rot = 90,              # rotate OTU labels for readability
  labsize = 5               # adjust label size if crowded
)

