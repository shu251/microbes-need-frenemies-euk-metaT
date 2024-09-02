# determine appropriate number of kmeans

# Load library
library(tidyverse)
library(broom)
# Load data
load("/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/kmeans-data-obj.RData", verbose = TRUE)


frenemies_26_clust <- kmeans(df_output, centers = 26)

summary(frenemies_26_clust)


broom::tidy(frenemies_26_clust)

augment_frenemies <- augment(frenemies_26_clust, df_output)
df_kmeans_frenemies <- (data.frame(augment_frenemies))

save(df_kmeans_frenemies, file = "/scratch/group/hu-lab/frenemies/euk-metaT-eukrhythmic-output/kmeansoutput-frenemies-metaT.RData")


cat("\n\nDone\n\n")
