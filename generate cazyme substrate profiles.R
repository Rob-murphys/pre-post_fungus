setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")

library(dplyr)
source("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/scripts/final scripts/dbcan4 output family cleaner.R")

### Metagenome profiles ####
#===========================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(), "/download/*dbsub.out", sep = ""))

# Reading in all data files
subs = lapply(filelist, read.table, header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = NULL)

substrate_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_substrates(subs[[i]])
  substrate_profiles = rbind(substrate_profiles, current_profile)
}

write.csv(substrate_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/substrate_profiles.csv")
saveRDS(substrate_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/substrate_profiles.RDS")
#### Metagenome profiles end ####


### Termitomyces profiles ####
#============================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(), "/termitomyces_download/*dbsub.out", sep = ""))

# Reading in all data files
subs = lapply(filelist, read.table, header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = NULL)

substrate_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_substrates(subs[[i]])
  substrate_profiles = rbind(substrate_profiles, current_profile)
}

write.csv(substrate_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_substrate_profiles.csv")
saveRDS(substrate_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_substrate_profiles.RDS")
#### Termitomyces profiles end ####



### M. natalensis profiles ####
#=============================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(), "/mNatalensis_download/*dbsub.out", sep = ""))

# Reading in all data files
subs = lapply(filelist, read.table, header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = NULL)

substrate_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_substrates(subs[[i]])
  substrate_profiles = rbind(substrate_profiles, current_profile)
}

write.csv(substrate_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_substrate_profiles.csv")
saveRDS(substrate_profiles, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_substrate_profiles.RDS")
#### M. natalensis profiles end ####