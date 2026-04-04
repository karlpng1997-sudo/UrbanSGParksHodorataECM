# ============================================================
# Script: ASV_to_OTU.R
# Purpose: Cluster ASVs into OTUs (97% similarity) using DECIPHER
# Input: ASV_table.csv (ASVs in columns, samples in rows)
# Output: OTU_table.csv, OTU_map.csv
# ============================================================

setwd("C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025")

# ---- 1. Load packages ----
if (!require("DECIPHER")) install.packages("DECIPHER", dependencies = TRUE)
if (!require("Biostrings")) install.packages("Biostrings", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)

library(DECIPHER)
library(Biostrings)
library(dplyr)

# ---- 2. Load ASV table ----
asv <- read.csv("ASV_table.csv", row.names = 1, check.names = FALSE)

cat("Loaded ASV table with", ncol(asv), "ASVs and", nrow(asv), "samples.\n")

# ---- 3. Extract ASV sequences ----
seqs <- DNAStringSet(colnames(asv))

# ---- 4. Compute distance matrix and cluster into OTUs ----
cat("Computing distance matrix (this may take a few minutes)...\n")
dna_dist <- DistanceMatrix(seqs, verbose = FALSE)

cat("Clustering sequences into OTUs at 97% similarity...\n")

otu_clusters <- DECIPHER::TreeLine(
  myDistMatrix=dna_dist,
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters")

# ---- 5. Create ASV-to-OTU map ----
otu_map <- data.frame(ASV = colnames(asv), OTU = otu_clusters$cluster)
write.csv(otu_map, "OTU_map.csv", row.names = FALSE)
cat("OTU map saved as OTU_map.csv\n")

# ---- 6. Collapse ASV counts by OTU ----
asv_t <- as.data.frame(t(asv))
asv_t$OTU <- otu_map$OTU

otu_table <- asv_t %>%
  group_by(OTU) %>%
  summarise(across(everything(), sum))

otu_table <- as.data.frame(otu_table)
rownames(otu_table) <- otu_table$OTU
otu_table$OTU <- NULL
otu_table <- t(otu_table)

write.csv(otu_table, "OTU_table.csv", row.names = TRUE)
cat("OTU abundance table saved as OTU_table.csv\n")

cat("✅ Done! You now have OTUs clustered at 97% similarity.\n")

rep_seq <- aggregate(ASV ~ OTU, otu_map, head, 1)
rep_seqs <- DNAStringSet(rep_seq$ASV)
names(rep_seqs) <- paste0("OTU_", rep_seq$OTU)

taxa <- assignTaxonomy(rep_seqs, "C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025/sh_general_release_dynamic_10.10.2017.fasta", multithread = TRUE)
write.csv(taxa, "taxa.csv", row.names = FALSE)






#Start from here, use the ps here to create new file




# ============================================================
# Script: phyloseq_merge_OTU_tax.R
# Purpose: Merge OTU abundance table and taxonomy table into a phyloseq object
# Input: OTU_table.csv, OTU_taxonomy.csv
# Output: phyloseq object (ps)
# ============================================================

setwd("C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025")

library(phyloseq)
library(tidyverse)

# ---- 2. Load OTU abundance table ----
# Rows = samples, Columns = OTUs
otu <- read.csv("OTU_table.csv", row.names = 1, check.names = FALSE)
otu_mat <- as.matrix(t(otu))  # transpose to make OTUs rows
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# ---- 3. Load taxonomy table ----
# Expect columns: OTU | Taxonomy (or one column with semicolon-separated taxonomy)
tax <- read.csv("taxa.csv", row.names = 1, check.names = FALSE)

# If taxonomy is in one column (e.g. "Kingdom;Phylum;Class;Order;Family;Genus;Species")
if (ncol(tax) == 1) {
  tax_split <- strsplit(tax[,1], ";")
  tax_mat <- do.call(rbind, lapply(tax_split, function(x) {
    length(x) <- 7  # pad missing levels
    return(x)
  }))
  colnames(tax_mat) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
} else {
  tax_mat <- as.matrix(tax)
}

# Match OTU IDs between tables
tax_mat <- tax_mat[rownames(tax_mat) %in% rownames(otu_mat), , drop = FALSE]
OTU <- prune_taxa(rownames(tax_mat), OTU)

# Create taxonomy table
TAX <- tax_table(as.matrix(tax_mat))

# ---- 4. Merge into phyloseq object ----
ps <- phyloseq(OTU, TAX)
ps

table(tax_table(ps)[, "Genus"])

# ---- Load sample metadata ----
meta <- read.csv("sample_data1.csv", row.names = 1, check.names = FALSE)

# Convert to phyloseq sample_data object
SAMP <- sample_data(meta)

# ---- Merge into phyloseq ----
ps <- merge_phyloseq(ps, SAMP)

# ---- Check result ----
ps
sample_data(ps)

plot_richness(ps, measures=c("Observed"), color='Location')
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
plot_bar(ps.prop, x = "Sample", y = "Abundance", fill ="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#######################################
# Summarize counts by Genus
ps_genus <- tax_glom(ps, taxrank = "Genus")

# Summarize counts by Species
ps_species <- tax_glom(ps, taxrank = "Species")

# Example for Genus level
otu_genus <- as.data.frame(otu_table(ps_genus))
tax_genus <- as.data.frame(tax_table(ps_genus))

# Add Genus names as a new column
otu_genus$Genus <- tax_genus$Genus

# Move Genus name to the first column
otu_genus <- otu_genus %>% relocate(Genus)

# View
head(otu_genus)

# Get taxonomy table
tax <- as.data.frame(tax_table(ps))

# Create combined label
tax$GenusSpecies <- paste(tax$Genus, tax$Species, sep = "_")

# Replace NAs and clean
tax$GenusSpecies <- gsub("_NA", "", tax$GenusSpecies)
tax$GenusSpecies <- gsub("NA_", "", tax$GenusSpecies)
tax$GenusSpecies <- paste("OTU", "_", seq_len(nrow(tax)), "_", tax$GenusSpecies)

# Assign as new row names for OTU table
otu <- as.data.frame(otu_table(ps))
rownames(otu) <- tax$GenusSpecies

otu <- t(otu)
# Save new genus-species table
write.csv(otu, "OTU_table_GenusSpecies.csv")

