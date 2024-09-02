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

# Import sample list
sample_merged <- read_delim(file = "input-docs/frenemies-pretximport.txt") 
# sample_merged
names(files) <- sample_merged$SAMPLE_REP  

cat("\nRunning check for all file names:\n")
names(files)

# Include taxonomic and functional IDs for the tx2gene step
cat("\nImporting tax and fxn annotation df\n")
tax_and_fxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, sep = "\t")

# TEST LINE 
# tax_and_fxn <- read.table("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, nrows = 250, sep = "\t")
# head(tax_and_fxn)

ptm <- proc.time()

tx2gene_in <- tax_and_fxn %>% 
  dplyr::mutate(SEQ_ID = stringr::str_remove(SequenceID, ".p[:digit:]$")) %>% 
  dplyr::select(SEQ_ID, TRANSCRIPT_ID = transcript_name)

# Run tximport step
# txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene_in)

txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene_in, txOut=TRUE, countsFromAbundance = "lengthScaledTPM")

rm(tx2gene_in)
rm(tax_and_fxn)

print(proc.time()-ptm)

cat("\nPrep txi and merged sample table for DEseq input\n")
# library(tidyverse)

# Import and align with the txi$counts output

# Convert to data frame.
sample_merged <- as.data.frame(sample_merged)
rownames(sample_merged) <- sample_merged$SAMPLE_REP
colnames(txi$counts) <- rownames(sample_merged)
colnames(txi$abundance) <- rownames(sample_merged)
colnames(txi$length) <- rownames(sample_merged)

cat("\nCheck colnames for txi are correct\n")
colnames(txi$counts)

cat("\nAnd they should match these rownames:\n")
rownames(sample_merged)


save(txi, sample_merged, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/tximport-sept-2024.RData")

cat("\nDONE\n")

# Saves txi object, as this step requires a lot of memory to run.
