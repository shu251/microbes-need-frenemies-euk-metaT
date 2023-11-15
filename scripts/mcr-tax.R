# Import MCR data, join with taxonomy information, summarize and save
library(tidyverse)
library(tximport)

load(file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-mcr.RData", verbose = TRUE)


mcr_counts_scaled <- makeCountsFromAbundance(
  as.matrix(txi_mcr_euk_annot$counts),
  as.matrix(txi_mcr_euk_annot$abundance),
  as.matrix(txi_mcr_euk_annot$length),
  countsFromAbundance = "scaledTPM")
# class(mcr_counts_scaled)
# head(mcr_counts_scaled)
mcr_df_scaled <- as.data.frame(mcr_counts_scaled)

names_orig <- colnames(mcr_df_scaled)
# names_orig
names_new_0 <- sub("_[^_]+$", "", names_orig)
# names_new_0
names_new_1 <- sub("MCR_\\d\\d_", "", names_new_0)
# names_new_1
colnames(mcr_df_scaled) <- names_new_1

mcr_mean_counts_df <- mcr_df_scaled %>%
  cbind(as.list(.) %>%
          Filter(is.numeric, .) %>% 
          split(names(.)) %>%
          lapply(as.data.frame) %>%
          lapply(rowMeans) %>%
          setNames(paste0("mean.", names(.)))) %>% 
  select(starts_with("mean"))

cat("\nImport tax and fxn information\n")
taxfxn <- read_delim(file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv")

mcr_annotated_mean <- dplyr::left_join(mcr_mean_counts_df %>%
                                         mutate(SequenceID = rownames(mcr_mean_counts_df)),
                                       taxfxn,
                                       by = "SequenceID")

mcr_df_for_tax <- mcr_annotated_mean %>% 
  select(full_classification, starts_with("mean"), SequenceID) %>% 
  distinct()

mcr_wtax_only <- mcr_df_for_tax %>% 
  select(-SequenceID) %>% 
  pivot_longer(cols = starts_with("mean"), names_to = "SAMPLE", values_to = "scaledTPM") %>% 
  group_by(SAMPLE, full_classification) %>% 
  summarise(SUM_scaledTPM = sum(scaledTPM)) 
 
save(mcr_mean_counts_df, mcr_wtax_only, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/mcr_df_tax_means.RData")
cat("\nDone\n")