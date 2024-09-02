# determine appropriate number of kmeans

# Load library
library(tidyverse)
library(broom)
# Load data
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/kmeans-data-obj.RData", verbose = TRUE)

# Modify input dataframe
#kmeans_input_pac <- mean_ctr_df %>% 
#  select(SequenceID, starts_with("mean.")) %>% 
#  distinct()

kclusts <-
  tibble(k = 5:35) %>%
  mutate(
    kclust = map(k, ~ kmeans(df_output, .x)),
    glanced = map(kclust, glance)
  )

output_kclusts <- kclusts %>%
  unnest(cols = c(glanced))

write_delim(output_kclusts, file = "output_kclusts.txt")
cat("\n\nSee output.\n\n")



cat("\n\nDone\n\n")
