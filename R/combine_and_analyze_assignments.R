# ------------------------------------------------------------------------------------------------------------------------#
# Compile and analyze kingdom-level assignments created by ./R/parse_kingodm_assignments.R
# This script loads, merges, and analyzes .csv files of kingdom-level taxonomic assignments of all studies
# Author: Geoffrey Zahn
# Requirements: tidyverse v 1.2.1
# ------------------------------------------------------------------------------------------------------------------------#

# packages ####
library(tidyverse)

# find .csv files ####
datapath <- "./data"
project_directories <- file.path(datapath,list.files(datapath))
assignment_files <- list.files(project_directories,pattern = "_kingdom.RDS",full.names = TRUE)

# load and combine kingdom-level data frames from each study ####
x=1
dfs=c()
for(i in assignment_files){
  assign(x=paste0(basename(i)),value=readRDS(i),envir = .GlobalEnv)
  dfs[x] <- paste0(basename(i))
  x=x+1
}


# join them all together with full_join() ####
obj_list = sapply(dfs,get)
full <- purrr::reduce(obj_list, full_join,by=names(which(table(unlist(sapply(obj_list,names))) == length(obj_list))),all=TRUE)


# Tidy the full data frame ####
names(full)[names(full)=="Var1"] <- "Kingdom"
names(full)[names(full)=="Freq"] <- "ESV_Count"
full$Kingdom <- unlist(purrr::map(str_split(full$Kingdom,"k__"),2))

full$DB[full$DB=="all_taxonomy"] <- "UNITE_Euk"
full$DB[full$DB=="fungi_taxonomy"] <- "UNITE"
full$DB[full$DB=="ncbi_taxonomy"] <- "UNITE+NCBI"


# Pick up here ......... ####

      # Look at more than just project and database; host, ecosystem, location, etc....

# summary stats

# plots

ggplot(full,aes(x=DB,y=Proportion,fill=Kingdom)) +
  geom_bar(stat="identity") +
  facet_wrap(~Main_Organism)


# models
fungi <- full %>% filter(Kingdom == "Fungi")
mod1 = aov(data=fungi, Proportion ~ DB + Project)
summary(mod1)
TukeyHSD(mod1)

