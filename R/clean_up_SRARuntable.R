library(tidyverse)

df = read_delim("data/SRA/Sra-Unfiltered-Runtable.txt",delim = "\t",guess_max = 10000)

# clean up lat_lon
df2 = df[!is.na(df$lat_lon),]
df3 = df2[df2$lat_lon != "not applicable",]
df3 = df3[grep(pattern="°",df3$lat_lon, invert = TRUE),]

# clean up assay type
df3 = df3[df3$Assay_Type == "AMPLICON",]

# remove unwanted columns
names(df3)
df3 <- df3 %>% select(-c("Center_Name","DATASTORE_filetype","DATASTORE_provider","DATASTORE_region",
                         "INSDC_center_alias","INSDC_center_name","INSDC_first_public","INSDC_last_update","INSDC_status","InsertSize"))

write_csv(df3,"./data/SRA/SraRunTable(2).csv",col_names = TRUE)

glimpse(df3)
table(df3$LibrarySelection)
