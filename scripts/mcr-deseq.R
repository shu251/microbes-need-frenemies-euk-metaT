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

tx2gene_in <- taxfxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$"))

all_transcripts <- as.character(tx2gene_in %>%
                                  select(transcript_name) %>% 
                                  .[["transcript_name"]])

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
cat("\n\nDESeq for Piccard vs. Von Damm\n\n\n")

txi_mcr_byfield <- subsetTxi(txi_mcr, vent_only_mcr, all_transcripts)

ds_tpm_field <- DESeqDataSetFromTximport(txi_mcr_byfield,
                                         colData = vent_only_mcr,
                                         design = ~0 + FIELDYR)

cat("\n\nUpreg in Von Damm is positive\n\n")
ds_tpm_field$FIELDYR <- factor(ds_tpm_field$FIELDYR, levels = c("Piccard2020", "VonDamm2020"))

ds_output_field <- DESeq2::DESeq(ds_tpm_field)
resultsNames(ds_output_field)
summary(ds_output_field)

save(txi_mcr_byfield, vent_only_mcr, ds_tpm_field, ds_output_field, file = "/scratch/group/hu-lab/frenemies/MCR_byfield_DESeq.RData")

##
cat("\n\nDESeq for Vent vs. non-vent\n\n\n")

txi_mcr_vent_nonvent <- subsetTxi(txi_mcr, mcr_no_tf, all_transcripts)

ds_tpm_bin_type <- DESeqDataSetFromTximport(txi_mcr_vent_nonvent,
                                            colData = mcr_no_tf,
                                            design = ~0 + TYPE_BIN)

cat("\n\nUpreg in vent is positive\n\n")
ds_tpm_bin_type$TYPE_BIN <- factor(ds_tpm_bin_type$TYPE_BIN, levels = c("Non-vent", "Vent"))

de_output_bin_type <- DESeq2::DESeq(ds_tpm_bin_type)
resultsNames(de_output_bin_type)
summary(de_output_bin_type)

save(txi_mcr_byfield, mcr_no_tf, ds_tpm_bin_type, de_output_bin_type, file = "/scratch/group/hu-lab/frenemies/MCR_bytype_DESeq.RData")

###
cat("\n\nDESeq for grazing vs. in situ\n\n\n")

txi_mcr_tf_insitu <- subsetTxi(txi_mcr, mcr_paired_tf, all_transcripts)

ds_tpm_exp <- DESeqDataSetFromTximport(txi_mcr_tf_insitu,
                                       colData = mcr_paired_tf,
                                       design = ~0 + EXP)

cat("\n\nUpreg in situ is positive\n\n")
ds_tpm_exp$EXP <- factor(ds_tpm_exp$EXP, levels = c("Tf", "insitu"))

ds_output_exp <- DESeq2::DESeq(ds_tpm_exp)
resultsNames(ds_output_exp)
summary(ds_output_exp)

save(txi_mcr_tf_insitu, mcr_paired_tf, ds_tpm_exp, ds_output_exp, file = "/scratch/group/hu-lab/frenemies/MCR_byexp_DESeq.RData")

## Compile
df_deseq_compiled <- (data.frame(ds_output_field) %>% 
                        add_column(DESEQ = "Von Damm vs. Piccard") %>% 
                        mutate(REGULATION = case_when(
                          log2FoldChange > 0 ~ "upregulated in Von Damm",
                          log2FoldChange < 0 ~ "upregulated in Piccard"
                        ))) %>% 
  bind_rows(data.frame(de_output_bin_type) %>% 
              add_column(DESEQ = "Vent vs. non-vent") %>% 
              mutate(REGULATION = case_when(
                log2FoldChange > 0 ~ "upregulated in vent",
                log2FoldChange < 0 ~ "upregulated in non-vent"
              ))) %>% 
  bind_rows(data.frame(ds_output_exp) %>% 
              add_column(DESEQ = "Grazing vs. in situ") %>% 
              mutate(REGULATION = case_when(
                log2FoldChange > 0 ~ "upregulated in situ",
                log2FoldChange < 0 ~ "upregulated in grazing Tf"
              ))) %>% 
  rownames_to_column(var = "SequenceID") %>% 
  SIGNIFICANT = case_when(
    pvalue <= 0.05 ~ "Significantly",
    TRUE ~ "Not significantly"
  )

stats_deseq <- df_deseq_compiled %>% 
  group_by(DESEQ, REGULATION, SIGNIFICANT) %>% 
  summarise(Total_num = n())

write_delim(stats_deseq, file = "output/stats_deseq.txt", delim = "\t")

save(df_deseq_compiled, stats_deseq, ds_output_exp, de_output_bin_type, ds_output_field, file ="/scratch/group/hu-lab/frenemies/DEseq_results_MCR.RData")

