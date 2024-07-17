library(dplyr)
library(effectsize)
library(vegan)
library(pairwiseAdonis)

#### Analysis of ECs 3.2.1.1 and 3.2.1.4 ####
#===========================================#

subs_rm = c("barley_beta_glucan", "chitosan", "dextran", "casein", "collagen") # substraights to remove

## Parsing AZCL metadata ##
azcl_meta = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/AZCL_data/azcl_enzymes.csv", stringsAsFactors = FALSE)

clean_azcl_meta = azcl_meta %>%filter(!AZCL_Substrate %in% subs_rm)

## Parsing data ##
azcl_data = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/AZCL_data/azcl_data.csv", stringsAsFactors = TRUE)


clean_azcl_data = azcl_data %>% filter(Caste1 != "soldier") %>% filter(!AZCL_Substrate %in% subs_rm) # removing specified substraights and soldier caste



clean_azcl_data_meta = clean_azcl_data %>%
  left_join(clean_azcl_meta[colnames(clean_azcl_meta) %in% c("AZCL_Substrate", "EC_number")], by = "AZCL_Substrate") %>%
  filter(EC_number %in% c("3.2.1.1", "3.2.1.4"))

#### Permanova of AZCL distance matrix ####
#=========================================#

azcl_data = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/AZCL_data/azcl_data.csv", stringsAsFactors = FALSE)

clean_azcl_data = azcl_data %>% filter(!(Caste1 == "soldier" |  AZCL_Substrate == "casein" | AZCL_Substrate == "collagen"))# removing all soldier rows and protease results
clean_azcl_data_wide = dcast(clean_azcl_data, Sample_ID ~ AZCL_Substrate, value.var = "Averages_technical_replicates")
rownames(clean_azcl_data_wide) = clean_azcl_data_wide$Sample_ID # setting sample names
clean_azcl_data_wide$Sample_ID = NULL # removing sample names column
clean_azcl_data_wide = clean_azcl_data_wide[!rownames(clean_azcl_data_wide) == 22,] # removing sample 22 as it only has one entry

# Make metadata df
azcl_metadata = clean_azcl_data[clean_azcl_data$AZCL_Substrate == "barley_beta_glucan",][,c(-2:-6, -11:-12)] # getting columns relevant to sample
azcl_metadata = azcl_metadata[order(azcl_metadata$Sample_ID),] # ordering dataframe by sample ID
azcl_metadata$Timepoint = as.factor(azcl_metadata$Timepoint)
azcl_metadata$Sample_ID = as.character(azcl_metadata$Sample_ID)

# generate relative abundances
rela_azcl = NULL
for(samp in unique(clean_azcl_data$Sample_ID)){
  df = clean_azcl_data %>% 
    filter(Sample_ID == samp) %>% 
    mutate(relative_abundance = Averages_technical_replicates/ sum(Averages_technical_replicates))
  rela_azcl = rbind(rela_azcl, df)
}

rela_azcl_wide = dcast(data = rela_azcl, Sample_ID~AZCL_Substrate, value.var = "relative_abundance")
rela_azcl_wide = rela_azcl_wide %>% column_to_rownames(var = "Sample_ID")


## Permanova ##

adonis.res = adonis2(clean_azcl_data_wide ~ Timepoint, data = azcl_metadata, permutations = 9999, method = "bray", by = "terms")

adonis_pairwise.res = pairwise.adonis2(clean_azcl_data_wide ~ Timepoint , data = azcl_metadata, method = "bray", nperm = 9999)
