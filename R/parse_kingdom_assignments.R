# ------------------------------------------------------------------------------------------------------------------------#
# Compare taxonomic assignments between databases for each study
# This script creates .csv files of kingdom-level taxonomic assignments for each taxonomic database within all studies
# Author: Geoffrey Zahn
# Requirements: phyloseq v 1.25.2; tidyverse v 1.2.1; vegan v 2.5-4
# ------------------------------------------------------------------------------------------------------------------------#

# packages ####
library(phyloseq)
library(tidyverse)
library(vegan)

# Find data ####
datapath <- "./data"
project_directories <- file.path(datapath,list.files(datapath))


# nested for-loop: iterates through each project, calculates raw and proportional numbers of assigned kingdoms for each taxonomic assignment method
# creates .csv file for each project

# For testing only:
# i=project_directories[1]
# j=list.files(file.path(i),pattern = "ps-UNITE",full.names = TRUE)[1]

for(i in project_directories){

  y=1
  for(j in list.files(file.path(i),pattern = "ps-UNITE",full.names = TRUE)){
    
    # load data ####
    ps <- readRDS(j)
    
    # make table of taxonomic assignments
    kingdom <- table(tax_table(ps)[,1])
    phylum <- table(tax_table(ps)[,2])
    
    # get name of database being used
    db_name = unlist(purrr::map(str_split(unlist(purrr::map(strsplit(j,"-"),2)),".RDS"),1))
    # get name of project being used
    proj_name = basename(i)
    
    # create data frame of kingdom
    df = as.data.frame(kingdom)
    df$DB <- db_name
    df$Project <- proj_name
    df$Proportion <- mutate(df,Freq / sum(df$Freq))
    df$Proportion <- df$Proportion$`Freq/sum(df$Freq)`
    
    # add alpha diversity measures ####
    all_diversity <- phyloseq::ntaxa(ps)
    ps_fungi <- subset_taxa(ps,Kingdom == "k__Fungi")
    fungal_diversity <- ntaxa(ps_fungi)
    
    df$Total_Richness <- all_diversity
    df$Fungal_Richness <- fungal_diversity
    
    ps_fungi_non_zero <- subset_samples(ps_fungi,(vegan::diversity(ps_fungi@otu_table) != 0))
    mean_shannon <- mean(diversity(ps_fungi@otu_table))
    all_shannon <- diversity(ps_fungi@otu_table)
    
    df$Fungal_Shannon_SD <- sd(all_shannon)
    df$Fungal_Shannon_Mean <- mean_shannon
    
    # add beta diversity measures ####
    dis <- vegdist(otu_table(ps_fungi))
    groups <- sample_data(ps_fungi)$BioProject
    mod <- betadisper(dis, groups)
    df$Beta_Dispersion_Mean <- mean(mod$distances)
    df$Beta_Dispersion_SD <- sd(mod$distances)
    
    
    
    
    # Add metadata such as location and environment
      # isolation source (most common from study)
    if("isolation_source" %in% names(sample_data(ps))){
    iso_source = table(sample_data(ps)$isolation_source)
    df$Main_Isolation_Source <- as.character(arrange(as.data.frame(iso_source),desc(Freq))$Var1[1])
    }
    
      # host (most common from study)
    if("Host" %in% names(sample_data(ps))){
    host = table(sample_data(ps)$Host)
    df$Main_Host <- as.character(arrange(as.data.frame(host),desc(Freq))$Var1[1])
    }
    
    # organism (most common from study)
    if("Organism" %in% names(sample_data(ps))){
      organism = table(sample_data(ps)$Organism)
      df$Main_Organism <- as.character(arrange(as.data.frame(organism),desc(Freq))$Var1[1])
    }
    
    
      # lat/lon (mean latitude and longitude from whole study)
      # N, E are positive ... S, W are negative
    
    latlon_list <- str_split(sample_data(ps)$lat_lon," ")
    
    South = purrr::map(latlon_list,2) == "S"
    South[South==TRUE] <- -1
    South[South==FALSE] <- 1
    lat <- as.numeric(purrr::map(latlon_list,1)) * South
    
    West = purrr::map(latlon_list,4) == "W"
    West[West==TRUE] <- -1
    West[West==FALSE] <- 1
    lon <- as.numeric(purrr::map(latlon_list,3)) * West
    
    df$Mean_Lat <- mean(lat)
    df$Mean_Lon <- mean(lon)  
    
    
    # assign data frame to object in .globalEnv
    assign(x = paste0("df_",proj_name,"_",y),value = df)
    
    
    df
  y=y+1
  }
  
  df1 <- get(paste0("df_",proj_name,"_",1))
  df2 <- get(paste0("df_",proj_name,"_",2))
  df3 <- get(paste0("df_",proj_name,"_",3))
  
  assign(x = paste0("df_",proj_name), value = rbind(df1,df2,df3))
  saveRDS(get(paste0("df_",proj_name)),file.path(i,paste0(proj_name,"_kingdom.RDS")))
}
