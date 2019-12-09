# ------------------------------------------------------------------------------------------------#
# Process ITS1-extracted reads and build a phyloseq object for a given SRA project
# Apply this script to each SRA Accession directory on FWD reads only
# Author: Geoffrey Zahn
# Requirements: dada2 v 1.9.0; phyloseq v 1.25.2; tidyverse v 1.2.1
# ------------------------------------------------------------------------------------------------#


# Load packages ####
library(phyloseq); packageVersion("phyloseq")
library(dada2); packageVersion("dada2")
library(tidyverse); packageVersion("tidyverse")


# define needed functions ####
getN <- function(x) sum(getUniques(x))


# set file paths ####
datapath <- "./data"
project_directories <- file.path(datapath,list.files(datapath))
ITS_directories <- file.path(project_directories,"seqs")

filtpath <- file.path(ITS_directories, "filtered") # Filtered files go into the filtered/ subdirectory

# make directory for filtered fqs if not already present
for(i in filtpath){
  if(!file_test("-d", i)) dir.create(i) 
}



# i=project_directories[1]
# Filter and trim inside for-loop ####

for(i in project_directories){
  
  # set up file paths and names
  fns <- sort(list.files(file.path(i,"seqs"), full.names = TRUE, pattern = "fastq.gz.ITS1.gz"))
  sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 1)

  filts <- file.path(i, "seqs/filtered", paste0(sample.names, "_filt.fastq.gz"))

  # filter and trim
  out <- filterAndTrim(fns, filts,
                       maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=16) 

  # learn error rates ####
  # Since some samples may have had zero reads pass QC, reassign filtFs and filtRs
  filts <- sort(list.files(file.path(i,"seqs/filtered"), full.names = TRUE))
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
  track <- cbind(out[,1], sapply(dada, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "nonchim")
  rownames(track) <- sample.names
  track = as.data.frame(track)
  track$filter.loss = (track[,1]-track[,2])/track[,1]
  write.csv(track, file = file.path(i,"Tracked_Filtration_Stats.csv"), row.names = TRUE)  
  
  # import metadata ####
  meta_path <- list.files(file.path(i),pattern = "SraRunTable",full.names = TRUE)
  meta = read.csv(meta_path,stringsAsFactors = FALSE)
  row.names(meta) <- meta$Run
  meta = meta[order(row.names(meta)),]  

  # remove missing samples excluded due to poor QC 
  row.names(seqtab.nochim) <- unlist(map(strsplit(row.names(seqtab.nochim), split = "_"),1))
  good.samples = (row.names(meta) %in%  row.names(seqtab.nochim))
  meta <- (meta[good.samples,])
  rm(good.samples)

  # Remove all seqs with fewer than 100 nucleotides
  keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
  seqtab.nochim <- seqtab.nochim[,keeper_esvs]
  
  # Assign taxonomy using each of the database options
  unite_fungi <- assignTaxonomy(seqtab.nochim, "./taxonomy/sh_general_release_dynamic_02.02.2019.fasta.gz", multithread=20) # make generic outgroups
  unite_all <- assignTaxonomy(seqtab.nochim, "./taxonomy/sh_general_release_dynamic_all_02.02.2019.fasta.gz", multithread=20) # make standard UNITE
  unite_ncbi <- assignTaxonomy(seqtab.nochim, "./taxonomy/sh_general_release_dynamic_plus_NCBI_02.02.2019.fasta.gz", multithread=20) # make standard UNITE
  
  # Hand off to Phyloseq and save RDS files ####
  ps.fungi <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                       sample_data(meta), 
                       tax_table(unite_fungi))
  saveRDS(ps.fungi, file = file.path(i,"ps_UNITE-fungi_taxonomy.RDS"))
  rm(list=c("ps.fungi","unite_fungi"))

  ps.all <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                          sample_data(meta), 
                          tax_table(unite_all))
  saveRDS(ps.all, file = file.path(i,"ps_UNITE-all_taxonomy.RDS"))
  rm(list=c("ps.all","unite_all"))
     
  ps.ncbi <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         sample_data(meta), 
                         tax_table(unite_ncbi))
  saveRDS(ps.ncbi,file = file.path(i,"ps_UNITE-ncbi_taxonomy.RDS"))
  rm(list=c("ps.ncbi","unite_ncbi"))
     
}

