# ------------------------------------------------------------------------------------------------------------------------#
# Analyses of kingdom-level assignments by project and database created by ./R/combine_and_analyze_assignments.R
# This script loads the full data frame to visualize and model taxonomic assignments
# Author: Geoffrey Zahn
# Requirements: tidyverse v 1.2.1
# ------------------------------------------------------------------------------------------------------------------------#

# Load packages ####
library(tidyverse)

# Load data ####
full = readRDS("./output/FULL_Kingdom_df.RDS")




# Summary Stats ####
names(full)

numcols <- c("ESV_Count","Proportion","Total_Richness","Fungal_Richness","Fungal_Shannon_SD","Fungal_Shannon_Mean","Beta_Dispersion_Mean",
             "Beta_Dispersion_SD","Mean_Lat","Mean_Lon")

for(i in numcols)(
  full[,i] <- as.numeric(full[,i])
)

summary(full$Proportion)



# Basic Plots ####

ggplot(full,aes(x=DB,y=Proportion,fill=Kingdom)) +
  geom_bar(stat="identity") +
  facet_wrap(~Main_Organism)


# models
names(full)
fungi <- full %>% filter(Kingdom == "Fungi")
mod1 = aov(data=fungi, Proportion ~ Project + DB)
mod2 = aov(data=fungi, Total_Richness ~ Project + DB)


summary(mod1)
summary(mod2)

TukeyHSD(mod1)










