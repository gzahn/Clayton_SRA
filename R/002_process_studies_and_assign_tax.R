library(tidyverse)
library(dada2)
library(phyloseq)

unite.ref <- "./taxonomy/sh_general_release_dynamic_16.10.2022.fasta.gz"
unite.euk.ref <- "./taxonomy/sh_general_release_dynamic_all_16.10.2022.fasta.gz"


ps_list_notax <- read_rds("./output/ps_list_notax.RDS")

# remove any faulty studies
ps_list_notax <- ps_list_notax[which(!is.na(ps_list_notax))]

# get accessions and label list elements
names(ps_list_notax) <- 
ps_list_notax %>% map(sample_data) %>% map("bio_project") %>% map(unique) %>% unlist()

full <- 
  ps_list_notax %>% 
  reduce(merge_phyloseq)

# UNITE Fungal
taxa_UNITE <- assignTaxonomy(full, unite.ref, multithread=16, tryRC = TRUE)
# UNITE + Euk
taxa_UNITE_Euk <- assignTaxonomy(full, unite.euk.ref, multithread=16, tryRC = TRUE)
