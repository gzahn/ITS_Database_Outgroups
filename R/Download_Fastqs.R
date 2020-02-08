# ---------------------------------------------------------------------------------------------------------#
# Download all fastq files associated with all SRA projects listed in "./metadata/Study_List.xlsx"
# Save only the forward reads, removing the reverse reads
# Author: Geoffrey Zahn
# Requirements: sra-toolkit v 2.9.2
# ---------------------------------------------------------------------------------------------------------#

starttime <- Sys.time()


# find directories and files ####
datapath <- "./data"
project_directories <- file.path(list.dirs(datapath,recursive = FALSE))


# For testing only  !!! Do NOT leave in final version !!!!
project_directories = project_directories[5]

runtables <- list.files(datapath,recursive = TRUE,pattern = "SraRunTable",full.names = TRUE) 

# For testing only  !!! Do NOT leave in final version !!!!
runtables = runtables[5]


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
  print(i)
  # Get frame of SRA Accessions and metadata
  accessions <- runtable$Run # remove limitation of 1:5 after test runs
  
  # nested loop iterates through Accession List within each project directory
  for(j in accessions){
    getfqs(SRR_ID = j)  
  }
  
  # remove reverse reads
  print("Removing Reverse Reads, if any")
  seqpath <- file.path(project_directories[x],"seqs")
  file.remove(list.files(seqpath,pattern = "_pass_2.fastq.gz",full.names = TRUE))
  
  x <- x+1
}

endtime <- Sys.time()
elapsedtime <- difftime(endtime, starttime,units = "hours")
print(elapsedtime)

