## Import salmon count files
# SLURM: metaT-tximport-gr-mcr.sm
# Must run scripts/create-samplelist.R to generate the sample_merged_txi.txt file first 

setwd("/home/skhu/microbes-need-frenemies-euk-metaT") # REMOTE:grace

# Load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(tximport)
library(readr)

num_threads <- getDTthreads()
# num_threads

# Get list of all salmon output files
files <- Sys.glob("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/salmon/*_quant/quant.sf")

files_gr_mcr <- files[!grepl("AXIAL", files)]
files_gr_mcr_insitu <- files_gr_mcr[!grepl("_Tf_", files_gr_mcr)]

# To get only MCR files
# mcr_files <- Sys.glob("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/salmon/MCR_*_quant/quant.sf")

# Import sample list
sample_merged <- read_delim(file = "/home/skhu/microbes-need-frenemies-euk-metaT/input-docs/frenemies-pretximport.txt") 
sample_merged_gr_mcr <- sample_merged %>% filter(SITE != "AXIAL" & EXP != "Tf")

cat("Check file names and sample names\n")
head(files_gr_mcr_insitu)
head(sample_merged_gr_mcr)

names(files_gr_mcr_insitu) <- sample_merged_gr_mcr$SAMPLE_REP  

# Include taxonomic and functional IDs for the tx2gene step
cat("\nImporting tax and fxn annotation df\n")
taxfxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

# TEST LINE 
# tax_and_fxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, nrows = 250, sep = "\t")
# head(tax_and_fxn)

ptm <- proc.time()

tx2gene_in <- taxfxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$")) %>% 
  dplyr::select(SEQ_ID, TRANSCRIPT_ID = transcript_name)

# Run tximport step
cat("Run tximport step")
txi_frenemies_mcr_gr <- tximport::tximport(files_gr_mcr_insitu, type = "salmon", tx2gene = tx2gene_in, txOut=TRUE, countsFromAbundance = "lengthScaledTPM")

print(proc.time()-ptm)

cat("\nPrep txi and merged sample table for DEseq input\n")

# Import and align with the txi$counts output
# Convert to data frame.
sample_merged_gr_mcr <- as.data.frame(sample_merged_gr_mcr)
rownames(sample_merged_gr_mcr) <- sample_merged_gr_mcr$SAMPLE_REP
colnames(txi_frenemies_mcr_gr$counts) <- rownames(sample_merged_gr_mcr)
colnames(txi_frenemies_mcr_gr$abundance) <- rownames(sample_merged_gr_mcr)
colnames(txi_frenemies_mcr_gr$length) <- rownames(sample_merged_gr_mcr)

cat("\nCheck colnames for txi are correct\n")
colnames(txi_frenemies_mcr_gr$counts)

cat("\nAnd they should match these rownames:\n")
rownames(sample_merged_gr_mcr)

cat("saving frenemies txi output..\n")
save(txi_frenemies_mcr_gr, sample_merged_gr_mcr, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-feb-2025.RData")

cat("\nDONE with first save of txi objects\n")
# Saves txi object, as this step requires a lot of memory to run.

cat("\nMake counts from:")
lenscaled_TPM <- makeCountsFromAbundance(
  as.matrix(txi_frenemies_mcr_gr$counts),
  as.matrix(txi_frenemies_mcr_gr$abundance),
  as.matrix(txi_frenemies_mcr_gr$length),
  countsFromAbundance = "lengthScaledTPM"
)
head(lenscaled_TPM)

# Go from matrix to data frame. Remove rows that only have zeroes
gr_mcr_lenscaled_TPM_df <- as.data.frame(lenscaled_TPM) %>% 
  filter(if_any(everything(.), ~. != 0))

glimpse(gr_mcr_lenscaled_TPM_df)

## Get mean across replicates
counts_df <- gr_mcr_lenscaled_TPM_df
names_orig <- colnames(counts_df)
names_new <- sub("_[^_]+$", "", names_orig)
colnames(counts_df) <- names_new

mean_gr_mcr_TPM_df <- counts_df %>% 
  cbind(as.list(.) %>%
          Filter(is.numeric, .) %>%
          split(names(.)) %>%
          lapply(as.data.frame) %>%
          lapply(rowMeans) %>%
          setNames(paste0("mean.", names(.)))) %>%
  select(starts_with("mean"))

save(lenscaled_TPM, gr_mcr_lenscaled_TPM_df, mean_gr_mcr_TPM_df, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/dfs_gr_mcr_feb2025.RData")

