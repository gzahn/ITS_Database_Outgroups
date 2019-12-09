# ------------------------------------------------------------------------------------------------#
# Process ITS1-extracted reads and build a phyloseq object for a given SRA project
# Apply this script to each SRA Accession directory on FWD reads only
# Author: Geoffrey Zahn
# Requirements: dada2 v 1.9.0; phyloseq v 1.25.2; tidyverse v 1.2.1
# ------------------------------------------------------------------------------------------------#


# Load packages ####
# library(DESeq2);packageVersion("DESeq2")
# library(DECIPHER)
library(phyloseq); packageVersion("phyloseq")
library(dada2); packageVersion("dada2")
library(tidyverse); packageVersion("tidyverse")


# File parsing ####

path <- "./ITS1" # rename path
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "fastq.gz.ITS1.gz"))
tail(fns)
sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 1)


# Filter and trim ####
filts <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fns, filts, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=16) 


no.pass = which(as.data.frame(out)$reads.out == 0)
passing_files = out[-no.pass,]


# learn error rates ####
# Since some samples may have had zero reads pass QC, reassign filtFs and filtRs
filts <- sort(list.files(filtpath, full.names = TRUE))
errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)

# Dereplication ####
derep <- derepFastq(filts, verbose=TRUE)

# Since some samples may have been removed (no reads passed QC), reassign sample.names
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)
names(derep) <- sample.names

# Sample inference ####
dadaFs <- dada(derep, err=errF, multithread=TRUE)

# Make a sequence table ####
seqtab <- makeSequenceTable(dadaFs)


# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# reassign "out" and "seqtab.nochim" to remove missing reads
out = out[as.data.frame(out)$reads.out > 0,]
seqtab.nochim <- seqtab.nochim[sample.names,]
dada = dadaFs[sample.names]

# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1], sapply(dada, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "./output/Sequence_Stats.csv", row.names = TRUE)  # rename path


# import metadata ####
meta = read.csv("./CoralFungiSGMetadata.csv") # rename to generic sra metadata

row.names(meta) <- meta$SampleID # make sure this is SRA compatible

# reorder metadata
meta = meta[order(row.names(meta)),]

# remove missing samples excluded due to poor QC 
row.names(seqtab.nochim) <- unlist(map(strsplit(row.names(seqtab.nochim), split = "_"),1))
good.samples = (row.names(meta) %in%  row.names(seqtab.nochim))
meta <- (meta[good.samples,])
rm(good.samples)


# Remove all seqs with fewer than 100 nucleotides
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# Assign Taxonomy ####

taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/sh_general_release_dynamic_s_01.12.2017_w_anthozoa.fasta", multithread=20) # make generic outgroups
taxa2 <- assignTaxonomy(seqtab.nochim, "./taxonomy/sh_general_release_dynamic_s_01.12.2017.fasta", multithread=20) # make standard UNITE


# Hand off to Phyloseq ####
ps.UNITE <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))

ps.OUTGROUP <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                        sample_data(meta), 
                        tax_table(taxa2))


# Save RDS objects for Phyloseq
saveRDS(ps.UNITE, file = "./output/ps_UNITE-only_taxonomy.RDS")  # rename path
saveRDS(ps.OUTGROUP, file = "./output/ps_UNITE-plus-outgroups_taxonomy.RDS") # rename path
