setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")

library(tibble)
library(dplyr)
library(ape)
library(stringr)
source("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/scripts/final scripts/dbcan4 output family cleaner.R")

#### Metagenome profiles ####
#===========================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(),"/download/*overview.txt", sep = ""))

# Reading in all data files
cazys = lapply(filelist, read.csv, header = TRUE, sep = "\t", stringsAsFactors=FALSE)

## cleaning and recovering only those rows with families identified using the logic of a custom function ##
family_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_familys(cazys[[i]])
  family_profiles = rbind(family_profiles, current_profile)
}

write.csv(family_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/family_profiles.csv")
saveRDS(family_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/family_profiles.RDS")
#### Metagenome profiles end ####

#### Termitomyces profiles ####
#=============================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(),"/termitomyces_download/*overview.txt", sep = ""))

# Reading in all data files
cazys = lapply(filelist, read.csv, header = TRUE, sep = "\t", stringsAsFactors=FALSE)

## cleaning and recovering only those rows with families identified using the logic of a custom function ##
family_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_familys(cazys[[i]])
  current_profile$sample_name = strsplit(basename(filelist[[i]]), "_overview")[[1]][1]
  family_profiles = rbind(family_profiles, current_profile)
}


write.csv(family_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_family_profiles.csv")
saveRDS(family_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_family_profiles.RDS")
#### Termitomyces profiles end ####


#### M. natalensis profiles ####
#=============================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(),"/mNatalensis_download/*overview.txt", sep = ""))

# Reading in all data files
cazys = lapply(filelist, read.csv, header = TRUE, sep = "\t", stringsAsFactors=FALSE)

## cleaning and recovering only those rows with families identified using the logic of a custom function ##
family_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_familys(cazys[[i]])
  current_profile$sample_name = strsplit(basename(filelist[[i]]), "_overview")[[1]][1]
  family_profiles = rbind(family_profiles, current_profile)
}

write.csv(family_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_family_profiles.csv")
saveRDS(family_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_family_profiles.RDS")
#### M. natalensis profiles end ####


## Simple family stats ##
#=======================#
count_df = as.data.frame(matrix(ncol = 2))
colnames(count_df) = c("family", "count")
count_df = count_df[-1,]
for(x in unique(family_profiles$family)){
  df = family_profiles %>%
    filter(family == x)
  count_df$family
  family = x
  count = nrow(df)/nrow(family_profiles)*100
  count_df = rbind(count_df, cbind(count, family))
}

## Summarized family stats ##

fam_groups = c("GH", "CBM", "PL", "CE", "GT", "AA")

for(x in fam_groups){
  group_df = count_df %>%
    filter(str_detect(family, x))
  print(c(x,sum(as.numeric(group_df$count))))
}

