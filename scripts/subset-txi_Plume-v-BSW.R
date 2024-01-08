setwd("/home/skhu/microbes-need-frenemies-euk-metaT")

library(DESeq2)
library(tidyverse)
library(tximport)
library(data.table)

load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-nov-2023.RData", verbose = TRUE)
# sample_merged_new <- sample_merged %>% 
# mutate(VENT = str_replace(VENT, "Mustard Stand", "MustardStand"))
# write_delim(sample_merged_new, file = "input-docs/sample_merged_txi.txt")

sample_merged <- read_delim("input-docs/sample_merged_txi.txt")

# Use same name ending with "_REP1" or "_REP2"

rm_plume <- c("Intl District plume bottom", "Intl District plume top", "Intl District plume +30m above top")

bsw_plume <- sample_merged %>% 
  filter(TYPE == "Plume" | TYPE == "Background") %>% 
  filter(!(SAMPLE_NAME %in% rm_plume)) %>% 
  filter(EXP == "insitu") %>% 
  select(sample = SAMPLE_REP)

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

## TESTING
# taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t", nrows = 300)
## TESING

tx2gene_in <- taxfxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$"))

euks_only <- as.character(tx2gene_in %>% 
                            filter(grepl("Eukaryota", full_classification)) %>% 
                            select(transcript_name) %>% 
                            .[["transcript_name"]])

euks_annot_only <- as.character(
  filter(tx2gene_in, grepl("Eukaryota", full_classification) & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(transcript_name) %>% 
    .[["transcript_name"]])

# save(all, all_no_tf, mcr_no_tf, mcr, mcr_paired_tf, axial, euks_only, euks_annot_only, file = "input-docs/objs-txi-subset.RData")
# cat("\n\nDONE\n\n")

# Subset txi directly
subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
{
  genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
  txi$abundance <- txi$abundance[genes, samples$sample]
  txi$counts <- txi$counts[genes, samples$sample]
  txi$length <- txi$length[genes, samples$sample]
  return(txi)
}

cat("\n\nStart txi subset\n\n")
# All samples, no grazing experiments
txi_euk_annot_plumebsw <- subsetTxi(txi, bsw_plume, euks_annot_only)
txi_euk_only_plumebsw <- subsetTxi(txi, bsw_plume, euks_only)

save(txi_euk_annot_plumebsw, txi_euk_only_plumebsw, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-plumebsw_samples.RData")

#
cat("\n\nCompleted txi subset, repeat for DESeq\n\n\n")
#
load(file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-plumebsw_samples.RData")
library(DESeq2)
library(tidyverse)
library(tximport)

# Import data
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-nov-2023.RData", verbose = TRUE)
load(file = "input-docs/objs-txi-subset.RData", verbose = TRUE)

# Import metadata
sample_merged <- read_delim("input-docs/sample_merged_txi.txt")

# Function to subset txi and compare vent vs. non-vent samples
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


de_bsw_plume_annotall <- deseq_bsw_plume(bsw_plume, euks_annot_only)
de_bsw_plume <- deseq_bsw_plume(bsw_plume, euks_only)

save(de_bsw_plume, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/background_vs_plume_DESeq.RData")
cat("\n\nDONE\n\n")
