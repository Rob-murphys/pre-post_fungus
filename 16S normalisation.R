library(vegan)
library(dplyr)
library(tibble)
library(stringr)
library(compositions)

data_handle = function(infile, outfile1){
  # Import data
  feature_table = as.data.frame(t(readRDS(infile)))

  # Data normalisation
  feature_table.norm = feature_table %>%
    sqrt(.) %>% # Squareoot standardization
    wisconsin(.) # Wisconson normalisation
  
  saveRDS(feature_table.norm, file = outfile1)
}

data_handle(infile = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/pre_post_feature_table_decontam.rds",
            outfile1 = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/pre_post_feature_table_normalised.rds")
