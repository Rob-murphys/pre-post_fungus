setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")

# Generating a filelist for all summary.txt files that are 3 subdirectories deep from the pwd to get EC numbers
filelist = Sys.glob(paste(getwd(), "/*/*overview.txt", sep = ""))

# Reading in all data files
ECs = lapply(filelist, read.csv, header = TRUE, sep = "\t", stringsAsFactors=FALSE)


# Correcting the EC number column and adding the sequence ID to each cazyme by splitting up the filepath and taking a parent directory as the seqID
fixed_EC = list()
for (i in 1:length(filelist)){
  df = ECs[[i]]
  df$EC. = sapply(df$EC., function(x){
    strsplit(as.character(x), "|", fixed = TRUE)[[1]][1]} ) # for column functions in each dataframe in the list split by "|" and return only the first part
  df$EC. = sapply(df$EC., function(x){
    strsplit(as.character(x), ":", fixed = TRUE)[[1]][1]} )
  df = df[!grepl("-", df$EC.),] # removing rows without EC numbers and those with incomplete ECs
  df$Gene.ID = sapply(df$Gene.ID, function(x){
    paste(strsplit(as.character(x), "_", fixed = TRUE)[[1]][1:3], collapse = "_")})
  fixed_EC = c(fixed_EC, list(df))
  fixed_EC[[i]]$seq_id = strsplit(filelist[i], "/", fixed = TRUE)[[1]][length(strsplit(filelist[i], "/", fixed = TRUE)[[1]])-1] # this final number controls which parent directory
}

# concatenating all the dataframes into one
concat_cazy = as.data.frame(do.call(rbind, fixed_EC)) # binding all cazy dataframes into one

# Writing this to file
write.csv(concat_cazy, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/EC_profiles.csv", row.names = FALSE)
saveRDS(concat_cazy, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/EC_profiles.csv")

# Writing each individual dataframe to a file
lapply(fixed_EC, function(x){
  write.csv(x, paste(as.character(unique(x$seq_id)), "_identified_cazymes.csv", sep = ""), row.names = FALSE)
})
