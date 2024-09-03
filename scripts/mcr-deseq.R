setwd("/home/skhu/microbes-need-frenemies-euk-metaT")

library(DESeq2)
library(tidyverse)
library(tximport)
library(data.table)

load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-MCR-sept-2024.RData", verbose = TRUE)
# txi_mcr and sample_merged_mcr

all_mcr <- sample_merged_mcr %>%
  select(sample = SAMPLE_REP)

vent_only_mcr <- sample_merged_mcr %>%
  filter(EXP == "insitu") %>% 
  filter(TYPE_BIN == "Vent") %>% 
  select(sample = SAMPLE_REP)

mcr_no_tf <- sample_merged_mcr %>% 
  filter(EXP == "insitu") %>% 
  select(sample = SAMPLE_REP)

# Of all the samples at MCR, which have paired in situ vs. Tf?
tmp_tf <- sample_merged_mcr %>% 
  filter(SITE == "MCR") %>% 
  # filter(VENT != "Background") %>% 
  # filter(VENT != "Plume") %>% 
  select(EXP, VENT) %>% 
  add_column(VAR = 1) %>% 
  pivot_wider(names_from = EXP, values_from = VAR, values_fn = sum) %>% 
  drop_na()

vent_wexperiments <- as.character(tmp_tf$VENT)

mcr_paired_tf <- sample_merged_mcr %>% 
  filter(VENT %in% vent_wexperiments) %>% 
  select(sample = SAMPLE_REP)

# mcr_paired_tf, mcr_no_tf, all_mcr, vent_only_mcr

taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")
# taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, nrows = 250, sep = "\t")
tx2gene_in <- taxfxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$"))
# head(tx2gene_in)
all_transcripts <- as.character(tx2gene_in %>%
                                  select(SEQ_ID) %>% 
                                  .[["SEQ_ID"]])

save(mcr_paired_tf, mcr_no_tf, all_mcr, vent_only_mcr, all_transcripts, file = "/scratch/group/hu-lab/frenemies/mcr-txi-subset-objects.RData")

# Function to subset TXI
subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
{
  genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
  txi$abundance <- txi$abundance[genes, samples$sample]
  txi$counts <- txi$counts[genes, samples$sample]
  txi$length <- txi$length[genes, samples$sample]
  return(txi)
}

##
cat("\n\n 1. DESeq for Piccard vs. Von Damm\n\n\n")

txi_mcr_output <- subsetTxi(txi_mcr, vent_only_mcr, all_transcripts)

subset_metadata_deseq <- sample_merged_mcr %>% 
  filter(SAMPLE_REP %in% (as.character(vent_only_mcr$sample))) # Change df here.

ds_tpm <- DESeqDataSetFromTximport(txi_mcr_output,
                                         colData = subset_metadata_deseq,
                                         design = ~0 + FIELDYR) # Change column header here

## Subsample: transcripts must appear in more than 1 sample and have a sum greater than 10
## "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts. This evolved out of use cases with Paul McMurdie's phyloseq package for metagenomic samples.

groupsize <- 1
keep <- rowSums(counts(ds_tpm) >= 10) >= groupsize
ds_tpm_filtered_0 <- ds_tpm[keep,]
ds_tpm_output_filtered_est <- estimateSizeFactors(ds_tpm_filtered_0, type = 'poscounts')

# Report
cat("\nStarted with ", dim(ds_tpm)[1], "observations. Filtering by 2 samples and 10 counts resulted in,", dim(ds_tpm_filtered_0)[1], ", which is", (100*(dim(ds_tpm_filtered_0)[1]/dim(ds_tpm)[1])), "% of the data.\n\n")

cat("\n\nUpreg in Von Damm is positive\n\n")
ds_tpm_output_filtered_est$FIELDYR <- factor(ds_tpm_output_filtered_est$FIELDYR, levels = c("Piccard2020", "VonDamm2020")) # Change factors

ds_output_result <- DESeq2::DESeq(ds_tpm_output_filtered_est)
resultsNames(ds_output_result)
summary(ds_output_result)

# Change naming schema
byfield <- (data.frame(results(ds_output_result)) %>% 
    add_column(DESEQ = "Von Damm vs. Piccard") %>% 
    mutate(REGULATION = case_when(
      log2FoldChange > 0 ~ "upregulated in Von Damm",
      log2FoldChange <= 0 ~ "upregulated in Piccard"
    )))
##
cat("\n\n 2. DESeq for Vent vs. non-vent\n\n\n")

txi_mcr_output <- subsetTxi(txi_mcr, mcr_no_tf, all_transcripts)
# 
subset_metadata_deseq <- sample_merged_mcr %>% 
  filter(SAMPLE_REP %in% (as.character(mcr_no_tf$sample))) # Change df here.

ds_tpm <- DESeqDataSetFromTximport(txi_mcr_output,
                                   colData = subset_metadata_deseq,
                                   design = ~0 + TYPE_BIN) # Change column header here
## Subsample: transcripts must appear in more than 1 sample and have a sum greater than 10
## "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts. This evolved out of use cases with Paul McMurdie's phyloseq package for metagenomic samples.
groupsize <- 1
keep <- rowSums(counts(ds_tpm) >= 10) >= groupsize
ds_tpm_filtered_0 <- ds_tpm[keep,]
ds_tpm_output_filtered_est <- estimateSizeFactors(ds_tpm_filtered_0, type = 'poscounts')

# Report
cat("\nStarted with ", dim(ds_tpm)[1], "observations. Filtering by 2 samples and 10 counts resulted in,", dim(ds_tpm_filtered_0)[1], ", which is", (100*(dim(ds_tpm_filtered_0)[1]/dim(ds_tpm)[1])), "% of the data.\n\n")

cat("\n\nUpreg in vent is positive\n\n")
ds_tpm_output_filtered_est$FIELDYR <- factor(ds_tpm_output_filtered_est$FIELDYR, levels = c("Non-vent", "Vent")) # Change factors

ds_output_result <- DESeq2::DESeq(ds_tpm_output_filtered_est)
resultsNames(ds_output_result)
summary(ds_output_result)

# Change naming schema
bybin_type <- (data.frame(results(ds_output_result)) %>% 
              add_column(DESEQ = "Vent vs. non-vent") %>% 
              mutate(REGULATION = case_when(
                log2FoldChange > 0 ~ "upregulated in vent",
                log2FoldChange <= 0 ~ "upregulated in non-vent"
              )))

###
cat("\n\n 3. DESeq for grazing vs. in situ\n\n\n")

txi_mcr_output <- subsetTxi(txi_mcr, mcr_paired_tf, all_transcripts)
 
subset_metadata_deseq <- sample_merged_mcr %>% 
  filter(SAMPLE_REP %in% (as.character(mcr_paired_tf$sample))) # Change df here.

ds_tpm <- DESeqDataSetFromTximport(txi_mcr_output,
                                   colData = subset_metadata_deseq,
                                   design = ~0 + EXP) # Change column header here
## Subsample: transcripts must appear in more than 1 sample and have a sum greater than 10
## "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts. This evolved out of use cases with Paul McMurdie's phyloseq package for metagenomic samples.
groupsize <- 1
keep <- rowSums(counts(ds_tpm) >= 10) >= groupsize
ds_tpm_filtered_0 <- ds_tpm[keep,]
ds_tpm_output_filtered_est <- estimateSizeFactors(ds_tpm_filtered_0, type = 'poscounts')

# Report
cat("\nStarted with ", dim(ds_tpm)[1], "observations. Filtering by 2 samples and 10 counts resulted in,", dim(ds_tpm_filtered_0)[1], ", which is", (100*(dim(ds_tpm_filtered_0)[1]/dim(ds_tpm)[1])), "% of the data.\n\n")

cat("\n\nUpreg in site is positive\n\n")
ds_tpm_output_filtered_est$FIELDYR <- factor(ds_tpm_output_filtered_est$FIELDYR, levels = c("Tf", "insitu")) # Change factors

ds_output_result <- DESeq2::DESeq(ds_tpm_output_filtered_est)
resultsNames(ds_output_result)
summary(ds_output_result)

# Change naming schema
by_exptype <- (data.frame(results(ds_output_result)) %>% 
                 add_column(DESEQ = "Grazing vs. in situ") %>% 
                 mutate(REGULATION = case_when(
                   log2FoldChange > 0 ~ "upregulated in situ",
                   log2FoldChange <= 0 ~ "upregulated in grazing Tf"
                 )))


## Compile
df_deseq_compiled <- byfield %>% 
  bind_rows(by_exptype) %>%
  bind_rows(bybin_type) %>%
  rownames_to_column(var = "SequenceID") %>% 
  mutate(SIGNIFICANT = case_when(
    pvalue <= 0.05 ~ "Significantly",
    TRUE ~ "Not significantly"
  ))

stats_deseq <- df_deseq_compiled %>% 
  group_by(DESEQ, REGULATION, SIGNIFICANT) %>% 
  summarise(Total_num = n())
# head(stats_deseq)

write_delim(stats_deseq, file = "stats_deseq.txt", delim = "\t")

save(df_deseq_compiled, stats_deseq, byfield, by_exptype, bybin_type, file ="/scratch/group/hu-lab/frenemies/DEseq_results_MCR.RData")

