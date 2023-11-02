# Get list of all salmon output files
files <- Sys.glob("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/salmon/*_quant/quant.sf")
# files

# Include metadata
sample_list <- read.csv("frenemies-metat-SAMPLELIST.csv")

# Create a sample list dataframe
sample_merged <- data.frame("Files"=files) %>%
  tidyr::separate(Files,sep="salmon/",into=c("STEM","Name")) %>%
  tidyr::separate(Name,sep="_quant",into=c("SampleName_rep","EXCESS")) %>%
  dplyr::select(SampleName_rep) %>%
  dplyr::mutate(Replicate = case_when(
    stringr::str_detect(SampleName_rep, "_1$") ~ "REP1",
    stringr::str_detect(SampleName_rep, "_2$") ~ "REP2"
  )) %>% 
  dplyr::mutate(SampleName = 
                  stringr::str_remove(SampleName_rep, "_[:digit:]$")) %>% 
  dplyr::left_join(sample_list, by = c("SampleName" = "SITE_NUM_FIELDYR_VENT_EXP_SAMPLEID")) %>% 
  tidyr::unite("Sample_rep", SampleName, Replicate, sep = "_", remove = FALSE)

names(files) <- sample_merged$Sample_rep  

# Use "sample_merged" again here
sample_merged_2 <- sample_merged %>% 
  mutate(TYPE_BIN = case_when(
    str_detect(VENT, "Plume|plume|IntlDist") ~ "Non-vent",
    str_detect(VENT, "BSW|Background|Transit") ~ "Non-vent",
    TRUE ~ "Vent")) %>% 
  mutate(TYPE = case_when(
    str_detect(VENT, "Plume|plume|IntlDist") ~ "Plume",
    str_detect(VENT, "BSW|Background|Transit") ~ "Background",
    TRUE ~ "Vent"))
# sample_merged_2
write_delim(sample_merged_2, file = "sample_merged_txi.txt")