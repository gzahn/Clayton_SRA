library(phyloseq)
library(vegan)
library(tidyverse)
library(easystats)
source("./R/plot_bar2.R")
source("./R/palettes.R")
theme_set(theme_bw())
unite <- readRDS("./Data/ps_UNITE.RDS")
euk <- readRDS("./Data/ps_UNITE_Euk.RDS")

# make sure order of taxa is same
if(identical(
  taxa_names(unite),
  taxa_names(euk)
)){
  unite_tax <- tax_table(unite) %>% as.data.frame()
  euk_tax <- tax_table(euk) %>% as.data.frame()  
}


# match up assignments
overall_agreement <- data.frame(
  kingdom = unite_tax$Kingdom == euk_tax$Kingdom,
  phylum = unite_tax$Phylum == euk_tax$Phylum,
  class = unite_tax$Class == euk_tax$Class,
  order = unite_tax$Order == euk_tax$Order,
  family = unite_tax$Family == euk_tax$Family,
  genus = unite_tax$Genus == euk_tax$Genus,
  species = unite_tax$Species == euk_tax$Species
)
overall_agreement[is.na(overall_agreement)] <- FALSE

head(overall_agreement)

# get metadata ####
ps_list_notax <- readRDS("./Data/ps_list_notax.RDS")

# remove any faulty studies
ps_list_notax <- ps_list_notax[which(!is.na(ps_list_notax))]

# get accessions and label list elements
names(ps_list_notax) <- 
  ps_list_notax %>% map(sample_data) %>% map("bio_project") %>% map(unique) %>% unlist()

# combine all studies into one ps object
full <- 
  ps_list_notax %>% 
  reduce(merge_phyloseq)

unite <- phyloseq(otu_table(full,taxa_are_rows = FALSE),
                  tax_table(unite),
                  sample_data(full))
euk <- phyloseq(otu_table(full,taxa_are_rows = FALSE),
                  tax_table(euk),
                  sample_data(full))



find_agreement <- function(x,y){
  agreement <- data.frame(
    kingdom = x$Kingdom == y$Kingdom,
    phylum = x$Phylum == y$Phylum,
    class = x$Class == y$Class,
    order = x$Order == y$Order,
    family = x$Family == y$Family,
    genus = x$Genus == y$Genus,
    species = x$Species == y$Species
  )
  agreement[is.na(agreement)] <- FALSE
  # agreement[,z] <- sam_data(x) %>% row.names()
  return(agreement)
}

# do as above, but for merged samples ####
# host
x <- unite %>% 
  merge_samples("host",fun = sum)
x <- subset_taxa(x,taxa_sums(x)>0)
x <- tax_table(x) %>% as.data.frame()
y <- euk %>% 
  merge_samples("host",fun = sum)
y <- subset_taxa(y,taxa_sums(y)>0)
y <- tax_table(y) %>% as.data.frame()

z <- merge_samples(unite,"host",fun=sum)
z@tax_table <- find_agreement(x,y) %>% tax_table()
colnames(z@tax_table) <- rank_names(unite)
host_ps <- z


# habitat
x <- unite %>% 
  merge_samples("habitat",fun = sum)
x <- subset_taxa(x,taxa_sums(x)>0)
x <- tax_table(x) %>% as.data.frame()
y <- euk %>% 
  merge_samples("habitat",fun = sum)
y <- subset_taxa(y,taxa_sums(y)>0)
y <- tax_table(y) %>% as.data.frame()

z <- merge_samples(unite,"habitat",fun=sum)
z@tax_table <- find_agreement(x,y) %>% tax_table()
colnames(z@tax_table) <- rank_names(unite)
habitat_ps <- z

# bioproject
x <- unite %>% 
  merge_samples("bio_project",fun = sum)
x <- subset_taxa(x,taxa_sums(x)>0)
x <- tax_table(x) %>% as.data.frame()
y <- euk %>% 
  merge_samples("bio_project",fun = sum)
y <- subset_taxa(y,taxa_sums(y)>0)
y <- tax_table(y) %>% as.data.frame()

z <- merge_samples(unite,"bio_project",fun=sum)
z@tax_table <- find_agreement(x,y) %>% tax_table()
colnames(z@tax_table) <- rank_names(unite)
bio_project_ps <- z


# examine agreement ####

# fix taxa_names
taxa_names(habitat_ps) <- colnames(otu_table(habitat_ps))
taxa_names(host_ps) <- colnames(otu_table(host_ps))
taxa_names(bio_project_ps) <- colnames(otu_table(bio_project_ps))


habitat_melt <- 
  habitat_ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt() %>% 
  mutate(Kingdom = Kingdom %>% as.logical(),
         Phylum = Phylum %>% as.logical(),
         Class = Class %>% as.logical(),
         Order = Order %>% as.logical(),
         Family = Family %>% as.logical(),
         Genus = Genus %>% as.logical(),
         Species = Species %>% as.logical()) %>% 
  pivot_longer(c(Kingdom,Phylum,Class,Order,Family,Genus,Species),
               names_to = "Rank",values_to = "Agreement")
host_melt <- 
  host_ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt() %>% 
  mutate(Kingdom = Kingdom %>% as.logical(),
         Phylum = Phylum %>% as.logical(),
         Class = Class %>% as.logical(),
         Order = Order %>% as.logical(),
         Family = Family %>% as.logical(),
         Genus = Genus %>% as.logical(),
         Species = Species %>% as.logical()) %>% 
  pivot_longer(c(Kingdom,Phylum,Class,Order,Family,Genus,Species),
               names_to = "Rank",values_to = "Agreement")

bio_project_melt <- 
  bio_project_ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt() %>% 
  mutate(Kingdom = Kingdom %>% as.logical(),
         Phylum = Phylum %>% as.logical(),
         Class = Class %>% as.logical(),
         Order = Order %>% as.logical(),
         Family = Family %>% as.logical(),
         Genus = Genus %>% as.logical(),
         Species = Species %>% as.logical()) %>% 
  pivot_longer(c(Kingdom,Phylum,Class,Order,Family,Genus,Species),
               names_to = "Rank",values_to = "Agreement")

habitat_melt %>% 
  group_by(Sample,Rank,Agreement) %>% 
  summarize(agreement = sum(Abundance)) %>% 
  mutate(Rank = factor(Rank,levels = rank_names(unite))) %>% 
  ggplot(aes(x=Rank,y=agreement,fill=Agreement)) +
  geom_col() +
  facet_wrap(~Sample) +
  scale_fill_manual(values = pal.discrete) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  labs(x="Taxonomic rank",y="Percent agreement")
ggsave("./output/manuscript_figs/agreement_by_habitat_barplot.png",dpi=400,width = 8,height = 6)  

host_melt %>% 
  group_by(Sample,Rank,Agreement) %>% 
  summarize(agreement = sum(Abundance)) %>% 
  mutate(Rank = factor(Rank,levels = rank_names(unite))) %>% 
  ggplot(aes(x=Rank,y=agreement,fill=Agreement)) +
  geom_col() +
  facet_wrap(~Sample) +
  scale_fill_manual(values = pal.discrete) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  labs(x="Taxonomic rank",y="Percent agreement")
ggsave("./output/manuscript_figs/agreement_by_host_barplot.png",dpi=400,width = 8,height = 6)  

# table version for manuscript details:
host_melt %>% 
  group_by(Sample,Rank,Agreement) %>% 
  summarize(agreement = sum(Abundance)) %>% 
  mutate(Rank = factor(Rank,levels = rank_names(unite))) %>% 
  filter(!Agreement) %>% 
  # mutate(agreement = 1-agreement) %>% 
  filter(Rank == "Kingdom") %>% 
  ungroup() %>% 
  select(Sample,agreement) %>% 
  rename("System"="Sample",
         "Proportion_Disagreement" = "agreement") %>% 
  arrange(Proportion_Disagreement) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic() # kingdom-level

bio_project_melt %>% 
  group_by(Sample,Rank,Agreement) %>% 
  summarize(agreement = sum(Abundance)) %>% 
  mutate(Rank = factor(Rank,levels = rank_names(unite))) %>% 
  ggplot(aes(x=Rank,y=agreement,fill=Agreement)) +
  geom_col() +
  facet_wrap(~Sample) +
  scale_fill_manual(values = pal.discrete) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  labs(x="Taxonomic rank",y="Percent agreement")
ggsave("./output/manuscript_figs/agreement_by_bioproject_barplot.png",dpi=400,width = 8,height = 6)  

# models ####
habitat_melt %>% 
  group_by(Sample,Rank,Agreement) %>% 
  summarize(agreement = sum(Abundance)) %>% 
  filter(Agreement) %>% 
  glm(data = .,
      formula = agreement ~ Sample * Rank) %>% 
  summary()

joined <- full_join(
  bio_project_melt,
  data.frame(
    Sample = unite@sam_data$bio_project,
    Host = unite@sam_data$host,
    Habitat = unite@sam_data$habitat
  )
)

metadata <- joined %>% 
  select(Sample,Host,Habitat) 
metadata <- 
metadata[seq(1,nrow(joined),by=nrow(kingdom_taxa)),] %>% unique.data.frame()
names(metadata)[1] <- "SRA Accession"
x <- metadata %>% 
  arrange(Host) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic()
saveRDS(x,"./output/manuscript_figs/metadata_table_kable.RDS")


kingdom_mod <- 
joined %>% 
  filter(Rank == "Kingdom",
         Agreement) %>% 
  glm(data=.,
      formula=Abundance ~ Host + Habitat)
summary(kingdom_mod)

kingdom_aov_mod <- 
  joined %>% 
  filter(Rank == "Kingdom",
         Agreement) %>% 
  aov(data=.,
      formula=Abundance ~ Host + Habitat)
summary(kingdom_aov_mod)
HSD <- kingdom_aov_mod %>% TukeyHSD(which = "Host")
saveRDS(kingdom_aov_mod,"./output/manuscript_figs/kingdom_aov_mod.RDS")

# HSD %>% plot()

# estimate_means(kingdom_mod,at = "Host") %>% plot()

# modelr::add_predictions(joined,kingdom_mod) %>% 
#   ggplot(aes(x=))

# which taxa were things actually being assigned when they disagreed? ####

kingdom_taxa <- data.frame(
  UNITE = unite@tax_table[,1] %>% str_remove_all("k__"),
  UNITE_Euk = euk@tax_table[,1] %>% str_remove_all("k__")
)
kingdom_taxa %>% 
  mutate(Euk_Kingdom = case_when(is.na(UNITE_Euk) ~ "unidentified",
                                 TRUE ~ UNITE_Euk),
         UNITE = case_when(is.na(UNITE) ~ "unidentified",
                           TRUE ~ UNITE)) %>% 
  mutate(disagree = UNITE != Euk_Kingdom) %>% 
  filter(disagree) %>% 
  select(UNITE,Euk_Kingdom) %>% table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>%  
  mutate(Euk_Kingdom = Euk_Kingdom %>% 
           str_remove("Eukaryota_kgd_") %>% 
           str_replace("unidentified","Unidentified")) %>% 
  mutate(Euk_Kingdom = factor(Euk_Kingdom,levels = Euk_Kingdom)) %>% 
  mutate(Freq = Freq/sum(Freq)) %>% 
  ggplot(aes(x=Euk_Kingdom,y=Freq)) +
  geom_col() +
  labs(x="'True' kingdom",y="Proportion of disagreements") +
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1))
ggsave("./output/manuscript_figs/kingdom-level_disagreement_taxonomy.png",dpi=400, height = 6, width = 8)



# clean up memory
rm(list=c("unite","unite_tax","joined","euk","euk_tax","ps_list_notax","x","y","z","full"))
