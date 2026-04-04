library(dada2)
library(dplyr)
path<- 'C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025'
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(151,151),
#                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(0, 0),    # no truncation for now
                     maxN = 0,
                     maxEE = c(5, 5),       # allow more expected errors
                     truncQ = 2,
                     rm.phix = FALSE,
                     compress = TRUE,
                     multithread = FALSE)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Make sequence table (ASV table)
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Check how many reads remain
sum(seqtab.nochim) / sum(seqtab)

# Write out representative sequences and table
write.csv(seqtab.nochim, "ASV_table.csv")

#Checks
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
sum(track[ ,"nonchim"])
sum(track[ ,"input"])
mean(track[ ,"nonchim"])

#Don't continue for OTUs
##############################################
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
seqtab <- (select(as.data.frame(seqtab), starts_with("CGTAACAAGGTTTCCGTAGG")))
seqtab <- do.call(cbind, lapply(seqtab, as.integer))
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025/sh_general_release_dynamic_10.10.2017.fasta", multithread=TRUE)
library(phyloseq)
library(Biostrings)
library(ggplot2)
x<-data.frame('Number'=1:20,'Place'=c("ZHNP Fungi", "ZHNP","ZHNP","ZHNP",
                                      "ZHNP","ZHNP","YP","YP","YP",
                                      "YP","YP","YP","SBG", "SBG Fungi",
                                      "SBG","SBG","SBG","Kallang",
                                      "Kallang Fungi","Kallang"))
rownames(x)<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(x), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
writeXStringSet(refseq(ps), "C:/Users/Karl/Desktop/karlfungi", format="fasta")

plot_richness(ps, measures=c("Observed"), color='Place')
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
p2=plot_ordination(ps.prop, ord.nmds.bray, color='Place',title="Bray NMDS")
p2+stat_ellipse() + theme_bw()
top30 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:30]#set all taxa found
ps.top30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top30 <- prune_taxa(top30, ps.top30)
plot_bar(ps, "Family", facet_grid=~Place)+ geom_bar(aes(fill=Family), stat="identity", position="stack")

setwd("C:/Users/Karl/Desktop")
library(ape)
library(phytools)
tree<-read.tree('phylotree.txt')
tree
tree[["tip.label"]]<-taxa_names(ps)
phy_tree(ps)<-tree
ps
plot_tree(ps, ladderize="left", label.tips="Family", color="Place")

