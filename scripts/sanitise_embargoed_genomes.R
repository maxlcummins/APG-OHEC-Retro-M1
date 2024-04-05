#Read in our packages
library(tidyverse)
library(xml2)

#Define our not in function
`%nin%` <- Negate(`%in%`)

#Sanitise inputs to remove embargoed genomes

#Read in the names of the quarantined samples
WW_quarantine <- read_csv("delims/purgelist.txt", col_names = 'name', show_col_types = FALSE)

#Generate list of potential gull quarantined samples (2019 cohort unpublished)
Gull_quarantine <- read_csv("delims/purgelist2.txt", col_names = 'name', show_col_types = FALSE)

#Combine lists
quarantine <- bind_rows(WW_quarantine, Gull_quarantine)

#Remove junk
rm(list = c("Gull_quarantine", "WW_quarantine"))

### TEA genome data
cleanme_qc_report_path <- "analysis/qc_report.txt"

cleanme_abritamr_path <- "analysis/abritamr_resistance.txt"

cleanme_mlst_path <- "analysis/mlst.txt"

cleanme_clermontyping_path <- "analysis/clermontyping.txt"

cleanme_ectyper_path <- "analysis/ectyper.txt"

cleanme_fimtyper_path <- "analysis/fimtyper.txt"

cleanme_IncF_RST_path <- "analysis/IncF_RST.txt"

cleanme_genotype_path <- "analysis/genotype.txt"

#Read in all our paths starting with TEA
for(i in ls(pattern = "cleanme_")){
        #Read in first df
        gene_df <- read_delim(file = get(i), delim = "\t", show_col_types = FALSE)
        #Create a prefix for colnames based on filename
        colname_prefix <- basename(get(i))
        #Trim prefix
        colname_prefix <- gsub("\\.txt", "", colname_prefix)
        if("Name" %in% colnames(gene_df)){
                gene_df <- gene_df %>% rename("name" = Name)
        }
        gene_df <- gene_df %>% mutate(name = str_replace(name, "\\.fasta$|\\.fna$|\\.fsa$|\\.fa$", ""))
        #Filter based on our list
        gene_df <- gene_df %>% filter(name %nin% quarantine$name)
        write_delim(x = gene_df, file = paste0("analysis/filtered/", colname_prefix, ".txt"), delim = "\t")
}

#Read in all our paths starting with TEA
for(i in list.files("analysis/filtered", full.names = T)){
        gene_df <- read_delim(file = i, delim = "\t", show_col_types = FALSE)
        print(i)
        print(head(gene_df))
}

#Sanitise metadata

#Read in metadata
metadata <- read_delim("delims/metadata.csv", 
                   delim = ",",
                   escape_double = FALSE,
                   trim_ws = TRUE)

#Filter QC report
metadata <- metadata %>% filter(name %nin% quarantine$name)

#Rewrite QC report
metadata %>% write_csv("delims/metadata.csv")
