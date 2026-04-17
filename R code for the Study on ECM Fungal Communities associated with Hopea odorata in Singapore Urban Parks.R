# R code for the "The ectomycorrhizal fungal community of Hopea odorata (Dipterocarpaceae) planted in tropical urban parks of Singapore"

# Step 1 ------ Extraction of Sequences ------

# Load relevant packages

library(dada2)
library(dplyr)

# Set Working Directory

path<- 'C:/Users/karlp/OneDrive/Desktop/Rcode'

# To ensure consistent results

set.seed(100)

# List files in path

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract Library ID Names

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(0, 0),    # no truncation for now
                     maxN = 0,
                     maxEE = c(2, 2),       
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = FALSE)
head(out)
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Make sequence table

seqtab <- makeSequenceTable(mergers)

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Check how many reads remain

sum(seqtab.nochim) / sum(seqtab)

# Write out representative sequences and table

write.csv(seqtab.nochim, "ASV_table.csv")

# Further checks on reads

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
sum(track[ ,"nonchim"])
sum(track[ ,"input"])
mean(track[ ,"nonchim"])

# Step 2 ------ Conversion of ASVs to OTUs ------

# Load relevant packages

if (!require("DECIPHER")) install.packages("DECIPHER", dependencies = TRUE)
if (!require("Biostrings")) install.packages("Biostrings", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)

library(DECIPHER)
library(Biostrings)

# Load ASV table

asv <- read.csv("ASV_table.csv", row.names = 1, check.names = FALSE)
cat("Loaded ASV table with", ncol(asv), "ASVs and", nrow(asv), "samples.\n")

# Extract ASV sequences

seqs <- DNAStringSet(colnames(asv))

# Compute distance matrix and cluster into OTUs

dna_dist <- DistanceMatrix(seqs, verbose = FALSE)

otu_clusters <- DECIPHER::TreeLine(
  myDistMatrix=dna_dist,
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters")

# Create ASV-to-OTU map

otu_map <- data.frame(ASV = colnames(asv), OTU = otu_clusters$cluster)
write.csv(otu_map, "OTU_map.csv", row.names = FALSE)

# Collapse ASV counts by OTU

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

rep_seq <- aggregate(ASV ~ OTU, otu_map, head, 1)
rep_seqs <- DNAStringSet(rep_seq$ASV)
names(rep_seqs) <- paste0("OTU_", rep_seq$OTU)

# Assign taxonomy to OTUs

taxa <- assignTaxonomy(rep_seqs, "C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025/sh_general_release_dynamic_10.10.2017.fasta", multithread = TRUE)
write.csv(taxa, "taxa.csv", row.names = FALSE)

# Load relevant packages

library(phyloseq)
library(tidyverse)

# Load OTU abundance table
# Rows = Samples, Columns = OTUs

# *To note: Add Sample Names assigned to Library ID Names in the OTU_table.csv and save as new file

otu <- read.csv("OTU_table.csv", row.names = 1, check.names = FALSE) 
otu_mat <- as.matrix(t(otu))  # transpose to make OTUs rows
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Load taxonomy table
# Expect columns: OTU | Taxonomy (or one column with semicolon-separated taxonomy)

# *To note: Add a left column with number IDs in ascending order to sorted OTUs in the taxa.csv and save as new file

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

# Merge OTU and TAX into phyloseq object

ps <- phyloseq(OTU, TAX)
ps

table(tax_table(ps)[, "Genus"])

# Load sample environmental data

meta <- read.csv("sample_environmental_data.csv", row.names = 1, check.names = FALSE)

# Convert to phyloseq sample_data object

SAMP <- sample_data(meta)

# Merge SAMP into phyloseq

ps <- merge_phyloseq(ps, SAMP)

# Check results

ps
sample_data(ps)

# Create plot to visualize observed alpha diversity

plot_richness(ps, measures=c("Observed"), color='Location')

# Create plot to visualize observed proportion of fungal diversity in samples

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
plot_bar(ps.prop, x = "Sample", y = "Abundance", fill ="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

# Provide new names for OTUs by Genus and Species

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

# Save new Genus-Species table

write.csv(otu, "OTU_table_GenusSpecies.csv")

# Step 3 ------ Creation of OTU barplots ------

# Load the CSV file as df

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

# Remove column name in df

colnames(df)[1] <- ""

# Save the renamed table

write.csv(df, "OTU_table_GenusSpecies_renamed.csv", row.names = FALSE, quote = FALSE)

# Load the table as new df

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

# Remove column name in df

colnames(df)[1] <- ""

# Save the renamed table

write.csv(df, "OTU_table_GenusSpecies_renumbered_by_genus_order.csv", row.names = FALSE, quote = FALSE)

# Load the CSV file with OTUs renumbered as df

df <- read.csv("OTU_table_GenusSpecies_renumbered_by_genus_order.csv", row.names = 1)

df_prop <- df / rowSums(df)

write.csv(df_prop, "OTU_table_GenusSpecies_proportions.csv", quote = FALSE)

# Load proportion or count table as df

df <- read.csv("OTU_table_GenusSpecies_proportions.csv", check.names = FALSE, row.names = 1)

# List identified ECM genera from FungalTraits

ecm_genera <- c("Russula", "Tomentella", "Lactarius", "Inocybe", "Clavulina", "Scleroderma")

# Identify ECM vs non-ECM columns

is_ecm <- sapply(colnames(df), function(x) any(sapply(ecm_genera, function(g) grepl(g, x, ignore.case = TRUE))))

# Create a new column summing all non-ECM columns

df$`Non-ECM fungi` <- rowSums(df[, !is_ecm, drop = FALSE])

# Keep only ECM columns + the new one

df_ecm_summary <- cbind(df[, is_ecm, drop = FALSE], df[, "Non-ECM fungi", drop = FALSE])

# Add "*" to "Non-ECM Fungi"

colnames(df_ecm_summary)[colnames(df_ecm_summary) == "Non-ECM fungi"] <- "*Non-ECM fungi"

# Save the result

write.csv(df_ecm_summary, "OTU_table_ECM_and_nonECM.csv", quote = FALSE)

# Load relevant packages

library(ggplot2)
library(tidyr)
library(dplyr)

# Read data
# *To note: In OTU_table_ECM_and_nonECM.csv, need to amend names of OTUs accordingly to assigned colours below

df <- read.csv("OTU_table_ECM_and_nonECM.csv", row.names = 1, check.names = FALSE)

# Convert row names (samples) to a column

df$Sample <- rownames(df)

# Convert from wide to long format

df_long <- pivot_longer(df, 
                        cols = -Sample, 
                        names_to = "OTU", 
                        values_to = "Proportion")

# Create barplot for ECM Fungal OTUs in each sample

group_colors <- c('Clavulina_OTU01' = 'red',
                  'Lactarius_OTU01' = 'palegreen1',
                  'Russula_OTU01' = 'pink1',
                  'Russula_OTU02' = 'palevioletred2',
                  'Russula_OTU03_violeipes' = 'palevioletred4',
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
                  'Tomentella_OTU13_parmastoana' = 'grey85',
                  'Inocybe_OTU01' = 'purple',
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

# Create barplot for ECM Fungal Genus in each sample

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

# Step 4 ------ Statistical Analysis ------

# Read CSV file

df <- read.csv("OTU_table_ECM_and_nonECM.csv")

# Replace all non-zero values with 1

df[df != 0] <- 1

# Save the modified table to calculate alpha diversity of ECM fungal OTUs in csv

write.csv(df, "OTU_table_ECM_and_nonECM_converted.csv", row.names = FALSE)

# Load relevant packages

library(lme4)        # glmer.nb
library(MuMIn)       # model selection (dredge, AICc, r.squaredGLMM)
library(DHARMa)      # diagnostics
library(performance) # extra checks (optional)
library(car)

# *To note: Load data with alpha diversity of ECM fungal OTUs measured from 'OTU_table_ECM_and_nonECM_converted.csv' and create new file with sample environmental data below

dat <- read.csv("ECMASV_GLMM_OTU.csv")

# Ensure expected columns exist

stopifnot(all(c("alpha","pH","clay","site") %in% names(dat)))

# Coerce types

dat$site <- as.factor(dat$site)

# Quick sanity checks

summary(dat[c("alpha","pH","clay","site")])

## Test Full model (NB-GLMM)

m_full <- glmer.nb(alpha ~ pH + clay + sand + (1 | site), data = dat)
cat("\n--- Full model summary ---\n")
print(summary(m_full))

##colinearity

model <- glmer.nb(alpha ~ pH + silt + clay + sand + (1 | site), data = dat)

# Compute VIFs

vif(model)

## Model selection (AICc, all subsets of fixed effects)
# Required for dredge

op_old <- options(na.action = "na.fail")
on.exit(options(op_old), add = TRUE)

ms <- dredge(m_full)                  # ranks by AICc
cat("\n--- Model selection table (AICc) ---\n")
print(ms)

best_model <- get.models(ms, 1)[[1]]  # best AICc model
cat("\n--- Best model formula ---\n")
print(formula(best_model))

## Best model summary (coefficients + p-values)

cat("\n--- Best model summary ---\n")
sbest <- summary(best_model)
print(sbest)

# Neat coefficient table (Estimate, SE, z, p)

cat("\n--- Coefficients (fixed effects) ---\n")
coef_tab <- as.data.frame(coef(summary(best_model)))
names(coef_tab) <- c("Estimate","Std.Error","z.value","Pr(>|z|)")
print(coef_tab)

## Goodness-of-fit: R2

cat("\n--- R² (marginal = fixed only; conditional = fixed + random) ---\n")
print(r.squaredGLMM(best_model))

## Diagnostics (DHARMa)

cat("\n--- DHARMa diagnostics ---\n")
sim_res <- simulateResiduals(best_model, n = 1000)

# Plots (opens a series of DHARMa diagnostic plots)

plot(sim_res)

# Uniformity (KS test), dispersion, zero-inflation, outliers

cat("\nUniformity test:\n");         print(testUniformity(sim_res))
cat("\nDispersion test:\n");         print(testDispersion(sim_res))
cat("\nZero-inflation test:\n");     print(testZeroInflation(sim_res))
cat("\nOutliers test:\n");           print(testOutliers(sim_res))

# Test CCA

# Load relevant packages

library(vegan)
library(tidyverse)

# *To note: Prepare data in R for CCA table using proportions of ECM Fungal OTUs from "OTU_table_ECM_and_nonECM.csv"

# Load OTU table

asv <- read.csv("OTU_table_ECM_and_nonECM_CCA.csv")

# Load environmental metadata

env <- read.csv("ECMASV_GLMM_OTU.csv")

# Merge by sample ID

merged <- inner_join(asv, env, by = "sample")

# Extract ASV matrix and env variables

asv_matrix <- merged %>% select(contains("OTU"))
env_data   <- merged %>% select(sample, pH, clay, sand, site)

# Remove samples with zero total abundance

keep <- rowSums(asv_matrix) > 0
asv_matrix <- asv_matrix[keep, ]
env_data   <- env_data[keep, ]

# Full partial CCA model

cca_full <- cca(asv_matrix ~ pH + clay + sand + Condition(site), data = env_data)

# Overall test

anova(cca_full)

# Test each term

anova(cca_full, by = "term")

# Test axes

anova(cca_full, by = "axis")

# Model selection (adjusted R²)

cca_null <- cca(asv_matrix ~ Condition(site), data = env_data)
cca_step <- ordiR2step(cca_null, scope = formula(cca_full),
                       direction = "both", data = env_data)

summary(cca_step)

# Plot selected model
ordiplot(cca_step, type = "n", main = "Partial CCA (selected model)")
points(cca_step, display = "sites", col = "blue", pch = 19)
text(cca_step, display = "species", cex = 0.7, col = "darkgreen")

# Fit env vectors (only pH, clay) and add significant ones
envfit_res <- envfit(cca_step, env_data[, c("pH", "clay")], permutations = 999)
plot(envfit_res, col = "red", p.max = 0.05)

# Comparison among variables

library(tidyverse)

# Load data

env <- read.csv("ECMASV_GLMM_OTU.csv")

# Calculate mean, sd, and se of alpha per site

alpha_stats <- env %>%
  group_by(site) %>%
  summarise(
    mean_alpha = mean(alpha, na.rm = TRUE),
    sd_alpha   = sd(alpha, na.rm = TRUE),
    n          = n(),
    se_alpha   = sd_alpha / sqrt(n)
  )

print(alpha_stats)

ggplot(alpha_stats, aes(x = site, y = mean_alpha, fill = site)) +
  geom_col() +
  labs(title = "Average Number of ECM Fungal OTUs per Location",
       x = "Park Location", y = "Number of ECM Fungal OTUs") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")   # removes redundant legend since site is on x-axis

library(car)       # Levene's test for homogeneity of variance
library(FSA)       # for Dunn's test (non-parametric post-hoc)

# Load data
env <- read.csv("ECMASV_GLMM_OTU.csv")

# Check assumptions

# Levene's test for homogeneity of variance
leveneTest(alpha ~ site, data = env)

# ANOVA (if assumptions ok)
anova_model <- aov(alpha ~ site, data = env)
summary(anova_model)

# Tukey post-hoc
TukeyHSD(anova_model)

# Kruskal-Wallis test (non-parametric)
kruskal.test(alpha ~ site, data = env)
kruskal.test(clay ~ site, data = env)
kruskal.test(sand ~ site, data = env)
kruskal.test(silt ~ site, data = env)
kruskal.test(pH ~ site, data = env)

# Load relevant packages

library(vegan)
library(tidyverse)

# Load data

asv <- read.csv("OTU_table_ECM_and_nonECM_CCA.csv") 

env <- read.csv("ECMASV_GLMM_OTU.csv")

# Merge datasets by sample ID

merged <- inner_join(asv, env, by = "sample")

# Extract OTU abundance data

asv_data <- merged %>%
  select(contains("OTU"))

# Ensure numeric and replace NA with 0

asv_data <- asv_data %>%
  mutate(across(everything(), ~ as.numeric(.))) %>%
  replace(is.na(.), 0)

# Remove samples with zero total abundance

row_sums <- rowSums(asv_data)
merged <- merged[row_sums > 0, ]
asv_data <- asv_data[row_sums > 0, ]

# Extract grouping factor

site <- merged$site

# ANOSIM

bray_dist <- vegdist(asv_data, method = "bray")

anosim_result <- anosim(bray_dist, site)

# Output

print(anosim_result)
plot(anosim_result)

# Step 5 ------ Creation of Bipartite Figures ------

library(RColorBrewer)
library(bipartite)
library(dplyr)
library(ggplot2)

# *To note: Load and prepare data for bipartite plots using "OTU_table_ECM_and_nonECM_converted.csv"

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

# Load relevant packages

library(readr)

# *To note: Load and prepare data for bipartite plots based on location using "OTU_table_ECM_and_nonECM_converted.csv"

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

# Save plot dimensions: width = 5000, height = 7000

# ------ The End ------