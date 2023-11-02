## Import salmon count files

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

# Include taxonomic and functional IDs for the tx2gene step

tax_and_fxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

# TEST LINE 
# tax_and_fxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, nrows = 250, sep = "\t")
head(tax_and_fxn)
ptm <- proc.time()

tx2gene_in <- tax_and_fxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$")) %>% 
  dplyr::select(SEQ_ID, TRANSCRIPT_ID = transcript_name)

# Run tximport step
txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene_in)

rm(tx2gene_in)
rm(tax_and_fxn)

print(proc.time()-ptm)


cat("\nPrep txi and merged sample table for DEseq input\n")
# library(tidyverse)

# Import and align with the txi$counts output
sample_merged_set <- read_delim("input-docs/sample_merged_txi.txt")
rownames(sample_merged_set) <- sample_merged_set$Sample_rep
rownames(sample_merged_set) <- colnames(txi$counts)

save(txi, sample_merged_set, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-oct-2023.RData")

# Saves txi object, as this step requires a lot of memory to run.