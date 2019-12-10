# ------------------------------------------------------------------------------------------------------------------------#
# Compare taxonomic assignments between databases for each study
# This script creates .csv files of kingdom-level taxonomic assignments for each taxonomic database within all studies
# Author: Geoffrey Zahn
# Requirements: phyloseq v 1.25.2; tidyverse v 1.2.1
# ------------------------------------------------------------------------------------------------------------------------#

# packages ####
library(phyloseq)
library(tidyverse)

# Find data ####
datapath <- "./data"
project_directories <- file.path(datapath,list.files(datapath))


# nested for-loop: iterates through each project, calculates raw and proportional numbers of assigned kingdoms for each taxonomic assignment method
# creates .csv file for each project
for(i in project_directories){

  y=1
  for(j in list.files(file.path(i),pattern = ".RDS",full.names = TRUE)){
    
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
  write.csv(get(paste0("df_",proj_name)),file.path(i,paste0(proj_name,"_kingdom.csv")),quote = FALSE,row.names = FALSE)
}
