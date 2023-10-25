setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling")

library(dplyr)
## Reading in the general data
filelist1 = Sys.glob(paste(getwd(), "/*_ncyc.csv", sep = ""))

# Reading in all data files
ncyc = lapply(filelist1, read.csv, header = TRUE, stringsAsFactors=FALSE)

ncyc_profile = do.call("rbind", ncyc)
ncyc_profile$Contig = sapply(ncyc_profile$Contig, function(x){
  y = strsplit(as.character(x), "_", fixed = TRUE)
  paste0(head(y[[1]],-1), collapse = "_")
})


write.csv(ncyc_profile, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile.csv")
saveRDS(ncyc_profile, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile-RDS")
