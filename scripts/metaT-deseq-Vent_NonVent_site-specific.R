setwd("/home/skhu/microbes-need-frenemies-euk-metaT/")
library(DESeq2)
library(tidyverse)
library(tximport)

# Import data
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-nov-2023.RData", verbose = TRUE)
load(file = "input-docs/objs-txi-subset.RData", verbose = TRUE)

# Import metadata
sample_merged <- read_delim("input-docs/sample_merged_txi.txt")
# unique(sample_merged$TYPE_BIN)

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
  ds_tpm_output_filtered <- estimateSizeFactors(ds_tpm_output_filtered_0, type = 'iterate')
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
# unique(sample_merged$TYPE_BIN)
## Apply to each sample subset

# From vnet field samples, what are the differences in vent vs non-vent?
cat("\nStart with MCR only\n\n")
cat("\n\nMCR samples, no grazing experiments:\n")
de_mcr_annot <- deseq_vent_novent(mcr_no_tf, euks_annot_only)

cat("\n\nMCR samples, no grazing experiments, with unannotated:\n")
de_mcr_unannot <- deseq_vent_novent(mcr_no_tf, euks_only)

cat("\n\nAxial samples:\n")
de_axial_annot <- deseq_vent_novent(axial, euks_annot_only)

save(de_mcr_annot, de_mcr_unannot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/mcr-insitu_vent_vs_novent_DESeq.RData")

save(de_axial_annot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/axial-insitu_vent_vs_novent_DESeq.RData")
