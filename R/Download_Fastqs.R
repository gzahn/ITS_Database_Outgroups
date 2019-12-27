# ---------------------------------------------------------------------------------------------------------#
# Download all fastq files associated with all SRA projects listed in "./metadata/Study_List.xlsx"
# Save only the forward reads, removing the reverse reads
# Author: Geoffrey Zahn
# Requirements: sra-toolkit v 2.9.2
# ---------------------------------------------------------------------------------------------------------#



# find directories and files ####
datapath <- "./data"
project_directories <- file.path(list.dirs(datapath,recursive = FALSE))

runtables <- list.files(datapath,recursive = TRUE,pattern = "SraRunTable",full.names = TRUE) 


# function for fastq-dump ####
getfqs <- function(SRR_ID){
  system(paste0("/usr/bin/sratoolkit.2.9.2-ubuntu64/bin//fastq-dump", # needs absolute filepath???
                " --outdir ",project_directories[x], "/seqs/",
                " --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ", 
                SRR_ID),
         wait = TRUE,intern = FALSE)
}


# use nested for-loops to do all projects iteratively ####


# main loop iterates through project directories
x <- 1
for(i in runtables){
  # assign run table
  runtable <- read.csv(i,stringsAsFactors = FALSE)
  # Get frame of SRA Accessions and metadata
  accessions <- runtable$Run # remove limitation of 1:5 after test runs
  
  # nested loop iterates through Accession List within each project directory
  for(j in accessions){
    getfqs(SRR_ID = j)  
  }
  
  # remove reverse reads
  seqpath <- file.path(project_directories[x],"seqs")
  file.remove(list.files(seqpath,pattern = "_pass_2.fastq.gz",full.names = TRUE))
  
  x <- x+1
}


