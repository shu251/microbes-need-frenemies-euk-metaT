setwd("/home/skhu/microbes-need-frenemies-euk-metaT")

library(tidyverse)
library(tximport)
library(data.table)

load(file = "/scratch/group/hu-lab/frenemies/dfs_mcr_mean-only_sept2024.RData")

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

# head(taxfxn)
# head(mean_mcr_TPM_df)

long_df <- mean_mcr_TPM_df %>%
  rowwise() %>%
  mutate(NUM_ZERO = sum(c_across(starts_with("mean.")) == 0)) %>%
  rownames_to_column(var = "SequenceID") %>%
  pivot_longer(cols = starts_with("mean."), values_to = "scaledTPM") %>%
  filter(scaledTPM > 0) %>% 
  separate(name, c("mean.field", "LIBRARY_NUM", "fieldyear", "LOCATION", "SAMPLETYPE", "SAMPLEID"), "_", 
           remove = FALSE) %>%
  select(-mean.field, -name) %>%
  mutate(VENT_FIELD = case_when(grepl("Piccard", fieldyear) ~ "Piccard",
                                grepl("VonDamm", fieldyear) ~ "Von Damm",
                                grepl("Axial", fieldyear) ~ "Axial",
                                grepl("Gorda", fieldyear) ~ "Gorda Ridge")) %>%
  mutate(VENT_BIN = case_when(
    (LOCATION == "Background" | LOCATION == "Plume" | LOCATION == "BSW") ~ "Non-vent",
    grepl("IntlDistrict", LOCATION) ~ "Non-vent",
    grepl("ASHES", LOCATION) ~ "Non-vent",
    TRUE ~ "Vent"
  ))

glimpse(long_df)

as_is <- c("Amoebozoa", "Apusozoa", "Excavata", "Hacrobia", "Archaeplastida")

long_df_annot <- long_df %>%
  left_join(taxfxn, by = "SequenceID") %>%
  separate(full_classification, c("Domain", "Supergroup", "Phylum", "Class", "Order", "Family", "Genus_spp"), sep = "; ", remove = FALSE) %>%
  mutate(SUPERGROUP_18S = case_when(
    Phylum == "Ciliophora" ~ "Alveolata-Ciliophora",
    Phylum == "Dinophyta" ~ "Alveolata-Dinoflagellata",
    # Phylum == "Perkinsea" ~ "Protalveolata",
    # Phylum == "Colponemidia" ~ "Protalveolata",
    # Phylum == "Chromerida" ~ "Protalveolata",
    Supergroup == "Alveolata" ~ "Other Alveolata",
    Supergroup %in% as_is ~ Supergroup,
    Supergroup == "Haptista" ~ "Hacrobia",
    Phylum == "Radiolaria" ~ "Rhizaria-Radiolaria",
    Phylum == "Cercozoa" ~ "Rhizaria-Cercozoa",
    (Supergroup == "Rhizaria" & Phylum != "Radiolaria" & Phylum != "Cercozoa") ~ "Rhizaria",
    Order == "Bigyra" ~ "Stramenopiles-Opalozoa;Sagenista",
    Class == "Ochromonadales" ~ "Stramenopiles-Ochrophyta",
    Supergroup == "Stramenopiles" ~ "Stramenopiles",
    Supergroup == "Opisthokonta" ~ "Opisthokonta",
    (is.na(Supergroup) | Supergroup == "Eukaryota incertae sedis") ~ "Unknown Eukaryota",
    TRUE ~ "Other-metaT only"))
# head(long_df_annot)
dim(long_df_annot);dim(long_df)

save(long_df_annot, long_df, file = "/scratch/group/hu-lab/frenemies/dfs_mcr_only.RData")
