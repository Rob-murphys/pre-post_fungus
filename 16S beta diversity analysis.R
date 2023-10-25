library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(pairwiseAdonis)

#### Import data ####
#===================#
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/metadata.csv")

metadata_noAncistro = metadata %>%
  filter(Species == "Macrotermes") %>%
  filter(BGI_seq_QC == "OK")

pre_post_feature_table = as.data.frame(readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/pre_post_feature_table_normalised.rds")) 

pre_post_feature_table_filt = pre_post_feature_table %>%
  arrange(as.numeric(row.names(.))) %>%
  filter(row.names(.) %in% metadata_noAncistro$Sample_ID)
#### Import data end ####


#### Bray curtis dissimilarity ####
#=================================#

bray_dis = vegdist(pre_post_feature_table_filt, method = "bray")
bray_dismat = as.matrix(pre_post_feature_table_filt) # as a matrix

#### Permanovas ####
#==================#

#Take output of pairwise adonis2 and makes it into a dataframe becuase staring at huge lists is not fun
adonis2df = function(adonis_res){
  res_df = NULL
  for(pair in 2:length(adonis_res)){
    new_row = cbind(from = strsplit(names(adonis_res[pair]), "_vs_")[[1]][1],
                    to = strsplit(names(adonis_res[pair]), "_vs_")[[1]][2],
                    R2 = adonis_res[pair][[1]][3][1,1], 
                    fval = adonis_res[pair][[1]][4][1,1], 
                    pval = adonis_res[pair][[1]][5][1,1])
    res_df = rbind(res_df, new_row)
  }
  return(res_df)
}

adonis_16S.res = adonis2(pre_post_feature_table_filt ~ Timepoint, by = "terms", data = metadata_noAncistro, method = "bray", permutations = 9999)

adonis_pairwise_16S.res = pairwise.adonis2(pre_post_feature_table_filt ~ Timepoint, data = metadata_noAncistro, method = "bray", nperm = 9999)
#### Bray curtis dissimilarity end ####
