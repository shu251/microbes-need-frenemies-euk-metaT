library(tidyverse)

# Import data
load(file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/all_samples_vent_metaT.RData", verbose = TRUE)

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

# Make long format data
cat("\n Make long \n")

long_df <- mean_counts_df %>% 
  rownames_to_column(var = "SequenceID") %>%
  rowwise() %>% 
  mutate(NUM_ZERO = sum(c_across(starts_with("mean.")) == 0)) %>% 
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


# Join with taxonomic and functional information
cat("\n\n Start joining and re-set taxa\n\n")

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

save(long_df, long_df_annot, file = "tmp_longdf.RData")


long_df_3000 <- long_df %>% sample_n(3000)
long_df_annot_3000 <- long_df_annot %>% sample_n(3000)

save(long_df_3000, long_df_annot_3000, file = "tmp_longdf_3000.RData")

cat("\nSave, done.")
