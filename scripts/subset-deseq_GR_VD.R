# Subset txi for Gorda Ridge and Von Damm samples:
setwd("/home/skhu/microbes-need-frenemies-euk-metaT")

library(DESeq2)
library(tidyverse)
library(tximport)
library(data.table)

# Import data
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-nov-2023.RData", verbose = TRUE)

sample_merged <- read_delim("input-docs/sample_merged_txi.txt")

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

# Subsample samplelist for TXI subset step

gr_vd <- sample_merged %>% 
  filter(SITE != "AXIAL") %>% 
  filter(EXP == "insitu") %>% 
  select(sample = SAMPLE_REP)

gr_vd_vent <- sample_merged %>% 
  filter(SITE != "AXIAL") %>% 
  filter(TYPE == "Vent") %>% 
  filter(EXP == "insitu") %>% 
  select(sample = SAMPLE_REP)

gr_vd_bsw_plume <- sample_merged %>% 
  filter(SITE != "AXIAL") %>% 
  filter(TYPE == "Plume" | TYPE == "Background") %>% 
  filter(EXP == "insitu") %>% 
  select(sample = SAMPLE_REP)


# Subsample list of genes:
tx2gene_in <- taxfxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$"))

euks_annot_only <- as.character(
  filter(tx2gene_in, grepl("Eukaryota", full_classification) & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(transcript_name) %>% 
    .[["transcript_name"]])

# Subset txi directly
# subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
# {
#   genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
#   txi$abundance <- txi$abundance[genes, samples$sample]
#   txi$counts <- txi$counts[genes, samples$sample]
#   txi$length <- txi$length[genes, samples$sample]
#   return(txi)
# }
# 
# 
# cat("\n\nStart txi subset\n\n")
# # All samples, no grazing experiments
# txi_gr_vd <- subsetTxi(txi, gr_vd, euks_annot_only)
# txi_gr_vd_vent <- subsetTxi(txi, gr_vd_vent, euks_annot_only)
# txi_gr_vd_bsw_plume <- subsetTxi(txi, gr_vd_bsw_plume, euks_annot_only)

# Completed Jan 9
# save(txi_gr_vd, txi_gr_vd_vent, txi_gr_vd_bsw_plume, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-gr_vd_all.RData")
# cat("\n\n Saved txi Robjects \n\n")


# Extract TXI values:
# txi <- txi_gr_vd
# library(tximport)
# 
# counts_scaled <- makeCountsFromAbundance(
#   as.matrix(txi$counts),
#   as.matrix(txi$abundance),
#   as.matrix(txi$length),
#   countsFromAbundance = "scaledTPM"
#   # countsFromAbundance = "lengthScaledTPM"
# )

# counts_df <- as.data.frame(counts_scaled)
# 
# names_orig <- colnames(counts_df)
# names_new <- sub("_[^_]+$", "", names_orig)
# colnames(counts_df) <- names_new
# 
# mean_counts_df <- counts_df %>%
#   cbind(as.list(.) %>%
#           Filter(is.numeric, .) %>%
#           split(names(.)) %>%
#           lapply(as.data.frame) %>%
#           lapply(rowMeans) %>%
#           setNames(paste0("mean.", names(.)))) %>% 
#   select(starts_with("mean"))

# Completed Jan 9
# save(mean_counts_df, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/gr_vd_samples_vent_metaT.RData")

# cat("\n\n Saved count df from txi\n\n")

cat("\nImport mean_counts_df\n")
load(file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/gr_vd_samples_vent_metaT.RData", verbose = TRUE)

# Process to long format
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


#
as_is <- c("Amoebozoa", "Apusozoa", "Excavata", "Hacrobia", "Archaeplastida")
#

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

save(long_df, long_df_annot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/gr_vd_longdf.RData")

cat("\n\nSaved GR and VD longformat data\n\n")


# Generate DE transcript information:
## Re import data for consistency.
cat("\nImport txi objects\n")

load(file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-gr_vd_all.RData", verbose = TRUE)

sample_merged <- read_delim("input-docs/sample_merged_txi.txt")

# Function to subset txi and generate DE
deseq_bsw_plume <- function(sample_set, gene_set){
  # First incorporate the txi subset fxn
  subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
  {
    genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
    txi$abundance <- txi$abundance[genes, samples$sample]
    txi$counts <- txi$counts[genes, samples$sample]
    txi$length <- txi$length[genes, samples$sample]
    return(txi)
  }
  # Run subset and sample merge re-set
  txi_output <- subsetTxi(txi, sample_set, gene_set)
  # This script will keep replacing tmp_sample_merged
  tmp_sample_merged <- sample_merged %>% 
    filter(SAMPLE_REP %in% as.character(sample_set$sample))
  rownames(tmp_sample_merged) <- tmp_sample_merged$SAMPLE_REP
  rownames(tmp_sample_merged) <- colnames(txi_output$counts)
  #
  # Import as DESeq object - use LIGHT in the design
  # DESeq
  ds_tpm_output <- DESeqDataSetFromTximport(txi_output,
                                            colData = tmp_sample_merged,
                                            design = ~0 + TYPE)
  # return(ds_tpm_output)
  # Further process DESeq
  groupsize <- 2 # Transcript to consider, must be in at least 3 samples
  keep <- rowSums(counts(ds_tpm_output) >= 10) >= groupsize # And have >= to 10 counts
  ds_tpm_output_filtered_0 <- ds_tpm_output[keep,]
  ds_tpm_output_filtered <- estimateSizeFactors(ds_tpm_output_filtered_0, type = 'poscounts')
  ###
  # Filtering stats:
  cat("\nStarted with ", dim(ds_tpm_output)[1], "observations. Filtering by 2 samples and 10 counts resulted in,", dim(ds_tpm_output_filtered)[1], ", which is", (100*(dim(ds_tpm_output_filtered)[1]/dim(ds_tpm_output)[1])), "% of the data.\n\n")
  ###
  #
  ## Positive log fold change == up regulated in Plume, compared to Background
  ds_tpm_output_filtered$TYPE <- factor(ds_tpm_output_filtered$TYPE, levels = c("Background", "Plume"))
  de_loc_output <- DESeq2::DESeq(ds_tpm_output_filtered)
  resultsNames(de_loc_output)
  summary(de_loc_output)
  cat("\n\nCompleted\n\n")
  return(de_loc_output)
}

de_bsw_plume <- deseq_bsw_plume(gr_vd_bsw_plume, euks_annot_only)

### Repeat above:
# Function to subset txi and compare vent vs. non-vent samples
deseq_vent_novent <- function(sample_set, gene_set){
  # First incorporate the txi subset fxn
  subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
  {
    genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
    txi$abundance <- txi$abundance[genes, samples$sample]
    txi$counts <- txi$counts[genes, samples$sample]
    txi$length <- txi$length[genes, samples$sample]
    return(txi)
  }
  # Run subset and sample merge re-set
  txi_output <- subsetTxi(txi, sample_set, gene_set)
  # This script will keep replacing tmp_sample_merged
  tmp_sample_merged <- sample_merged %>% 
    filter(SAMPLE_REP %in% as.character(sample_set$sample))
  rownames(tmp_sample_merged) <- tmp_sample_merged$SAMPLE_REP
  rownames(tmp_sample_merged) <- colnames(txi_output$counts)
  #
  # Import as DESeq object - use LIGHT in the design
  # DESeq
  ds_tpm_output <- DESeqDataSetFromTximport(txi_output,
                                            colData = tmp_sample_merged,
                                            design = ~0 + TYPE_BIN)
  # return(ds_tpm_output)
  # Further process DESeq
  groupsize <- 2 # Transcript to consider, must be in at least 3 samples
  keep <- rowSums(counts(ds_tpm_output) >= 10) >= groupsize # And have >= to 10 counts
  ds_tpm_output_filtered_0 <- ds_tpm_output[keep,]
  ds_tpm_output_filtered <- estimateSizeFactors(ds_tpm_output_filtered_0, type = 'poscounts')
  ###
  # Filtering stats:
  cat("\nStarted with ", dim(ds_tpm_output)[1], "observations. Filtering by 2 samples and 10 counts resulted in,", dim(ds_tpm_output_filtered)[1], ", which is", (100*(dim(ds_tpm_output_filtered)[1]/dim(ds_tpm_output)[1])), "% of the data.\n\n")
  ###
  #
  ## Positive log fold change == up regulated in Vent, compared to Non-vent
  ds_tpm_output_filtered$TYPE_BIN <- factor(ds_tpm_output_filtered$TYPE_BIN, levels = c("Non-vent", "Vent"))
  de_loc_output <- DESeq2::DESeq(ds_tpm_output_filtered)
  resultsNames(de_loc_output)
  summary(de_loc_output)
  cat("\n\nCompleted\n\n")
  return(de_loc_output)
}

de_gr_vd <- deseq_vent_novent(gr_vd, euks_annot_only)

## Generated both Differential expression lists.
# de_bsw_plume
# de_gr_vd

# Process DE transcript signatures
library(DESeq2)
library(tidyverse)

## GR and Von Damm vent vs. non-vent:
results_all_fields <- DESeq2::results(de_gr_vd, alpha = 0.05)

plot_vent_v_nonvent <- data.frame(results_all_fields) %>% 
  mutate(REGULATION = case_when(
    log2FoldChange > 0 ~ "upregulated in vent",
    log2FoldChange < 0 ~ "upregulated in non-vent"
  ),
  SIGNIFICANT = case_when(
    pvalue <= 0.05 ~ "Significantly",
    TRUE ~ "Not significantly"
  )) %>% 
  ggplot(aes(x = baseMean, y = log2FoldChange, color = SIGNIFICANT)) +
  geom_point(stat = "identity") +
  scale_x_log10() +
  theme_classic() +
  scale_color_manual(values = c("#878787", "#d73027")) +
  labs(title = mcols(results_all_fields)$description[2])

allfields_vent_v_nonvent_transcripts <- data.frame(results_all_fields) %>% 
  mutate(REGULATION = case_when(
    log2FoldChange > 0 ~ "upregulated in vent",
    log2FoldChange < 0 ~ "upregulated in non-vent"
  ),
  SIGNIFICANT = case_when(
    pvalue <= 0.05 ~ "Significantly",
    TRUE ~ "Not significantly"
  )) %>% 
  # filter(SIGNIFICANT == "Significantly") %>% 
  rownames_to_column(var = "SequenceID")

cat("\n\nTable listing the number of transcripts in each category for GR and VD vent vs. non-vent\n")
allfields_vent_v_nonvent_transcripts %>% 
  group_by(REGULATION, SIGNIFICANT) %>% 
  summarise(Number_of_transcripts = n())


# Repeat for bsw and plume:
results_bsw_plume <- DESeq2::results(de_bsw_plume, alpha = 0.05)

plot_bsw_v_plume <- data.frame(results_bsw_plume) %>% 
  mutate(REGULATION = case_when(
    log2FoldChange > 0 ~ "upregulated in plume",
    log2FoldChange < 0 ~ "upregulated in background"
  ),
  SIGNIFICANT = case_when(
    pvalue <= 0.05 ~ "Significantly",
    TRUE ~ "Not significantly"
  )) %>% 
  ggplot(aes(x = baseMean, y = log2FoldChange, color = SIGNIFICANT)) +
  geom_point(stat = "identity") +
  scale_x_log10() +
  theme_classic() +
  scale_color_manual(values = c("#878787", "#d73027")) +
  labs(title = mcols(results_bsw_plume)$description[2])

plume_bsw_transcripts <- data.frame(results_bsw_plume) %>% 
  mutate(REGULATION = case_when(
    log2FoldChange > 0 ~ "upregulated in plume",
    log2FoldChange < 0 ~ "upregulated in background"
  ),
  SIGNIFICANT = case_when(
    pvalue <= 0.05 ~ "Significantly",
    TRUE ~ "Not significantly"
  )) %>% 
  # filter(SIGNIFICANT == "Significantly") %>% 
  rownames_to_column(var = "SequenceID")
cat("\n\nTable listing the number of transcripts in each category for GR and VD plume vs. background\n")
plume_bsw_transcripts %>% 
  group_by(REGULATION, SIGNIFICANT) %>% 
  summarise(Number_of_transcripts = n())

save(plot_vent_v_nonvent, allfields_vent_v_nonvent_transcripts, plot_bsw_v_plume, plume_bsw_transcripts, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/DESeq-output-GR_VD.RData")