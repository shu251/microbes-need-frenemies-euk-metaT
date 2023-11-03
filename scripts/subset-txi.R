setwd("/home/skhu/microbes-need-frenemies-euk-metaT")

library(DESeq2)
library(tidyverse)
library(tximport)
library(data.table)
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-oct-2023.RData", verbose = TRUE)

all <- sample_merged_set %>%
  select(sample = Sample_rep)

mcr <- sample_merged_set %>% 
  filter(SITE == "MCR") %>% 
  select(sample = Sample_rep)

axial <- sample_merged_set %>% 
  filter(SITE == "AXIAL") %>% 
  select(sample = Sample_rep)

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

tx2gene_in <- taxfxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$"))
# head(tx2gene_in)
euks_only <- as.character(tx2gene_in %>% 
                            filter(grepl("Eukaryota", full_classification)) %>% 
                            select(transcript_name) %>% 
                            .[["transcript_name"]])

euks_annot_only <- as.character(
  filter(tx2gene_in, grepl("Eukaryota", full_classification) & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(transcript_name) %>% 
    .[["transcript_name"]])

save(all, mcr, axial, euks_only, euks_annot_only, file = "input-docs/objs-txi-subset.RData")

# Subset txi directly
subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
{
  genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
  txi$abundance <- txi$abundance[genes, samples$sample]
  txi$counts <- txi$counts[genes, samples$sample]
  txi$length <- txi$length[genes, samples$sample]
  return(txi)
}

txi_euk_annot <- subsetTxi(txi, all, euks_annot_only)

txi_euk_only <- subsetTxi(txi, all, euks_only)

txi_mcr_euk_annot <- subsetTxi(txi, mcr, euks_annot_only)

txi_axial_euk_annot <- subsetTxi(txi, axial, euks_annot_only)

save(txi_euk_annot, txi_euk_only, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-allsamples.RData")
save(txi_mcr_euk_annot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-mcr.RData")
save(txi_axial_euk_annot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-axial.RData")


# txi_ca <- subsetTxi(txi, ca_only, euks_only)
# 
# tmp_sample_merged <- sample_merged %>% 
#   filter(Sample_rep %in% as.character(ca_only$sample))
# rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
# rownames(tmp_sample_merged) <- colnames(txi_ca$counts)
# 
# # Compare euphotic vs. subeuphotic in coastal California
# ## Includes Port of LA and Catalina
# ds_tpm_ca_light <- DESeqDataSetFromTximport(txi_ca,
#                                             colData = tmp_sample_merged,
#                                             design = ~0 + LIGHT)
# 
# save(ds_tpm_ca_light, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/ca-deseq.RData")