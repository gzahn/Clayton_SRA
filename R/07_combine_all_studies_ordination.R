# packages and setup ####
library(tidyverse)
library(phyloseq)
library(patchwork)
library(microbiome)
library(skimr)
source("./R/plot_bar2.R")
source("./R/palettes.R")

theme_set(theme_minimal())

# collect file paths 
path <- "./output/output_ps_data"
UNITE_list <- list.files(path,pattern = "UNITE_ps_object.RDS",full.names = TRUE)
EUK_list <- list.files(path,pattern = "UNITE_Euk_ps_object.RDS",full.names = TRUE)
StudyIDs <- UNITE_list %>% str_remove_all("./output/output_ps_data/") %>% str_remove_all("_UNITE_ps_object.RDS")

UNITES <- map(UNITE_list,readRDS)
EUKS <- map(EUK_list,readRDS)

full_UNITE_ps <- reduce(UNITES,merge_phyloseq)
full_EUK_ps <- reduce(EUKS,merge_phyloseq)

getmeta <- function(x){sample_data(x) %>% meta()}
full_meta <- map(UNITES,getmeta)
map(full_meta,glimpse)


# ordinate
nmds <- full_EUK_ps %>% 
  ordinate(method = "NMDS")

ord_by_study <- plot_ordination(physeq = full_EUK_ps, ordination = nmds,color="BioProject") +
  coord_cartesian(xlim = c(-75,-25),ylim = c(-75,-20)) +
  scale_color_viridis_d(end=.8)

ord_by_feature <- plot_ordination(physeq = full_EUK_ps, ordination = nmds,color="Organism") +
  coord_cartesian(xlim = c(-75,-25),ylim = c(-75,-20))

ord_by_study / ord_by_feature


# Add "DatabaseMethod" and merge everything

full_UNITE_ps@sam_data$Database <- "UNITE"
full_EUK_ps@sam_data$Database <- "UNITE+Euk"
sample_names(full_EUK_ps) <- paste0(sample_names(full_EUK_ps),"_UNITE_Euk")
everything <- merge_phyloseq(full_UNITE_ps,full_EUK_ps)



# Ordinate
full_nmds <- everything %>% 
  # tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS")

full_ord_plot <- plot_ordination(everything, full_nmds, color="BioProject")
full_ord_plot
full_ord_plot +
  coord_cartesian(xlim = c(-12,0),ylim = c(-8,2)) +
  facet_wrap(~Database)
