setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")

library(dplyr)
library(reshape2)

## Loading in abundance data ##
abundance_df = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/anvio_coverage_stats/abundance_anvio.csv")

#### Family data ####
#===================#

## Loading in cazyme data ##
family_profiles = readRDS("family_profiles.RDS")

## merging the two ##
family_profiles_abun = family_profiles %>%
  left_join(abundance_df, by = c("contig" = "contig"))



## Removing "_indexed" from the sample name ##
family_profiles_abun = family_profiles_abun %>% 
  mutate(sample_name = str_remove(sample_name, "_indexed"))

## Writing as a long dataframe

write.csv(family_profiles_abun, "family_profiles_abun.csv", row.names = FALSE)
saveRDS(family_profiles_abun, "family_profiles_abun.rds")

## wriring as wide dataframe ##

family_profiles_abun_wide = dcast(data = family_profiles_abun, sample_name ~ family, value.var = "rel_abun", fun.aggregate = sum)

write.csv(family_profiles_abun_wide, "family_profiles_abun_wide.csv", row.names = FALSE)
saveRDS(family_profiles_abun_wide, "family_profiles_abun_wide.rds")

#### Family data  end ####




#### Substrate data ####
#======================#

## Loading in cazyme data ##
substrate_profiles = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/substrate_profiles.RDS")

## merging the two ##
substrate_profiles_abun = substrate_profiles %>%
  left_join(abundance_df, by = c("Gene ID" = "contig"))



## Removing "_indexed" from the sample name ##
substrate_profiles_abun = substrate_profiles_abun %>% 
  mutate(sample_name = str_remove(sample_name, "_indexed"))

## Writing as a long dataframe

write.csv(substrate_profiles_abun, "substrate_profiles_abun.csv", row.names = FALSE)
saveRDS(substrate_profiles_abun, "substrate_profiles_abun.rds")

## wriring as wide dataframe ##

substrate_profiles_abun_wide = dcast(data = substrate_profiles_abun, sample_name ~ Substrate, value.var = "rel_abun", fun.aggregate = sum)

write.csv(substrate_profiles_abun_wide, "substrate_profiles_abun_wide.csv", row.names = FALSE)
saveRDS(substrate_profiles_abun_wide, "substrate_profiles_abun_wide.rds")

#### Substrate data  end ####



#### Nitrogen cycling ####
#========================#

ncyc_profile = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile-RDS")

ncyc_profile_abun = ncyc_profile %>%
  left_join(abundance_df, by = c("Contig" = "contig"))

## Removing "_indexed" from the sample name ##
ncyc_profile_abun = ncyc_profile_abun %>% 
  mutate(sample_name = str_remove(sample_name, "_indexed"))

## Writing as a long dataframe

write.csv(ncyc_profile_abun, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun.csv", row.names = FALSE)
saveRDS(ncyc_profile_abun, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun.rds")

## wriring as wide dataframe ##

ncyc_profile_abun_wide = dcast(data = ncyc_profile_abun, sample_name ~ gene_name, value.var = "rel_abun", fun.aggregate = sum)

write.csv(ncyc_profile_abun_wide, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun_wide.csv", row.names = FALSE)
saveRDS(ncyc_profile_abun_wide, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun_wide.rds")

#### Nitrogen cycling end ####