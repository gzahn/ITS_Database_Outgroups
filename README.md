# ITS_Database_Outgroups
Comparing "standard" UNITE taxonomic assignments with UNITE+Outgroups in published fungal ITS1 studies


**Workflow:**
1. Download sequence data from all SRA accessions listed in "./metadata/Study_List.xlsx" keeping only FWD reads
    + "./R/Download_Fastqs.R"
2. Extract ITS1 region
    + "./R/Extract_ITS1.R"
3. Process each project with DADA2, assign taxonomy both ways, save phyloseq objects
    + "./R/make_phyloseq_objects.R"



All analyses conducted with:
+ R verion 3.4.4
+ R-studio version 
+ GNOME Terminal 3.18.3 using VTE version 0.42.5 +GNUTLS
