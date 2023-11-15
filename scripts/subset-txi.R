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
all <- sample_merged %>%
  select(sample = SAMPLE_REP)

all_no_tf <- sample_merged %>% 
  filter(EXP == "insitu") %>% 
  select(sample = SAMPLE_REP)

# Of all the samples at MCR, which have paired in situ vs. Tf?
tmp_tf <- sample_merged %>% 
  filter(SITE == "MCR") %>% 
  # filter(VENT != "Background") %>% 
  # filter(VENT != "Plume") %>% 
  select(EXP, VENT) %>% 
  add_column(VAR = 1) %>% 
  pivot_wider(names_from = EXP, values_from = VAR, values_fn = sum) %>% 
  drop_na()

vent_wexperiments <- as.character(tmp_tf$VENT)

mcr <- sample_merged %>% 
  filter(SITE == "MCR") %>% 
  select(sample = SAMPLE_REP)

mcr_no_tf <- sample_merged %>% 
  filter(SITE == "MCR") %>% 
  select(sample = SAMPLE_REP)

mcr_paired_tf <- sample_merged %>% 
  filter(SITE == "MCR") %>% 
  filter(VENT %in% vent_wexperiments) %>% 
  select(sample = SAMPLE_REP)

axial <- sample_merged %>% 
  filter(SITE == "AXIAL") %>% 
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

save(all, all_no_tf, mcr_no_tf, mcr, mcr_paired_tf, axial, euks_only, euks_annot_only, file = "input-docs/objs-txi-subset.RData")
cat("\n\nDONE\n\n")

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
txi_euk_annot <- subsetTxi(txi, all_no_tf, euks_annot_only)
txi_euk_only <- subsetTxi(txi, all_no_tf, euks_only)

txi_mcr_euk_annot <- subsetTxi(txi, mcr_no_tf, euks_annot_only)
txi_mcr_euk_annot_paired_exps <- subsetTxi(txi, mcr_paired_tf, euks_annot_only)

txi_axial_euk_annot <- subsetTxi(txi, axial, euks_annot_only)

save(txi_euk_annot, txi_euk_only, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-allsamples.RData")

save(txi_mcr_euk_annot, txi_mcr_euk_annot_paired_exps, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-mcr.RData")

save(txi_axial_euk_annot, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/txi-axial.RData")

