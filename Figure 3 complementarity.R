setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")

library(tibble)
library(dplyr)
library(reshape2)
library(circlize)
#================================#
## Read in and formate the data ##
#================================#

## Metadata ##
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/sample_metadata.csv", stringsAsFactors = FALSE)
metadata$Matriline = as.factor(metadata$Matriline)
metadata$Timepoint = as.factor(metadata$Timepoint)
metadata = metadata[order(metadata$Sample),] #ordering alphabetically


## Metagenomes ##
family_profiles_abun = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/family_profiles_abun.rds") %>%
  mutate(present = ifelse(rel_abun > 0, 1, 0)) %>%
  left_join(metadata, by = c("sample_name" = "Sample")) %>%
  select(-c("Gene.ID", "EC.", "HMMER", "dbCAN_sub", "DIAMOND", "X.ofTools", "contig", "mean_cov"))

  
## Termitomyces ##
termitomyces_family_profiles = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_family_profiles.RDS") %>%
  mutate(present = 1)

termitomyces_family_profiles_wide = dcast(termitomyces_family_profiles, family ~ sample_name, value.var = "present", fun.aggregate = sum)

termitomyces_family_profiles_wide = as.data.frame(t(termitomyces_family_profiles_wide %>%
  column_to_rownames(var = "family") %>%
  mutate(present = 1))) %>%
  slice(5)
inFungus_cazyme = colnames(termitomyces_family_profiles_wide)

#wide_fungus_cazy = wide_fungus_cazy %>% dplyr::slice(x = 5)

## M. natalensis ##

termite_family_profile = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_family_profiles.RDS")


termiteCazymeSummary = unique(termite_family_profile$family)
termite_family_profile$present = 1
termite_family_profile_wide = as.data.frame(t(dcast(termite_family_profile, family ~ sample_name, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "family")))

inTermite_cazyme = colnames(termite_family_profile_wide)



#### Pre timepoint plot ####
#==========================#

# Adding ECs not present in the gut microbiome be ARE present in the termite and fungus to the dataframe and setting them to 0
preFamily_profiles_abun = family_profiles_abun %>%
  dcast( Timepoint ~ family, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "Timepoint")

preFamily_profiles_abun[,setdiff(inTermite_cazyme,colnames(preFamily_profiles_abun))] = 0
preFamily_profiles_abun[,setdiff(inFungus_cazyme,colnames(preFamily_profiles_abun))] = 0
preFamily_profiles_abun[preFamily_profiles_abun > 0] = 1 # setting those present in the micoribome to 1


# Adding ECs not present in the termite  be ARE present in the gut microbiome and fungus to the dataframe and setting them to 0
termite_family_profile_wide[,setdiff(colnames(preFamily_profiles_abun), inTermite_cazyme)] = 0
termite_family_profile_wide[,setdiff(inFungus_cazyme, inTermite_cazyme)] = 0
termite_family_profile_wide[termite_family_profile_wide > 0] = 2 # setting those present in the termite to 2

# Adding ECs not present in the fungus  be ARE present in the gut microbiome and termite to the dataframe and setting them to 0
termitomyces_family_profiles_wide[,setdiff(colnames(preFamily_profiles_abun), inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[,setdiff(inTermite_cazyme, inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[termitomyces_family_profiles_wide > 0] = 3 # setting those present in the fungus to 3


comb_data = rbind(preFamily_profiles_abun, termite_family_profile_wide, termitomyces_family_profiles_wide)

comb_data_reorder = as.data.frame(t(comb_data)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present)) %>%
  select(!c("POST", "FIELD"))


cols1 = colorRamp2(c(0.5, 1, 2, 3), c("white", "#55C8FA", "#8FD485", "white"))

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 3/complementarity_cazyme_pre.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder, col = cols1, cluster = T, track.height = 0.5)
dev.off()
circos.clear()

#### Pre timepoint plot end ####


#### Post timepoint plot ####
#===========================#

post_family_profiles_abun = family_profiles_abun %>%
  #mutate(present = ifelse(Timepoint != "POST", 0, present)) %>% 
  dcast( Timepoint ~ family, value.var = "present", fun.aggregate = sum) %>%
  #filter(!Timepoint %in% c("PRE", "FIELD")) %>%
  column_to_rownames(var = "Timepoint")

post_family_profiles_abun[,setdiff(inTermite_cazyme,colnames(post_family_profiles_abun))] = 0
post_family_profiles_abun[,setdiff(inFungus_cazyme,colnames(post_family_profiles_abun))] = 0
post_family_profiles_abun[post_family_profiles_abun > 0] = 1

termite_family_profile_wide[,setdiff(colnames(post_family_profiles_abun), inTermite_cazyme)] = 0
termite_family_profile_wide[,setdiff(inFungus_cazyme, inTermite_cazyme)] = 0
termite_family_profile_wide[termite_family_profile_wide > 0] = 2

termitomyces_family_profiles_wide[,setdiff(colnames(post_family_profiles_abun), inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[,setdiff(inTermite_cazyme, inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[termitomyces_family_profiles_wide > 0] = 3


comb_data_post = rbind(post_family_profiles_abun, termite_family_profile_wide, termitomyces_family_profiles_wide)

comb_data_reorder = as.data.frame(t(comb_data_post)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present)) %>%
  select(!c("PRE", "FIELD"))

cols1 = colorRamp2(c(0.5, 1, 2, 3), c("white", "#227EEE", "#8FD485", "#F3D84F"))

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 3/complementarity_cazyme_post.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder, col = cols1, cluster = T, track.height = 0.5)
dev.off()
circos.clear()
#### Post timepoint plot end ####


#### Field timepoint plot ####
#============================#

field_family_profiles_abun = family_profiles_abun %>%
  #mutate(present = ifelse(Timepoint != "FIELD", 0, present)) %>% 
  dcast( Timepoint ~ family, value.var = "present", fun.aggregate = sum) %>%
  #filter(!Timepoint %in% c("PRE", "POST")) %>%
  column_to_rownames(var = "Timepoint")

field_family_profiles_abun[,setdiff(inTermite_cazyme,colnames(field_family_profiles_abun))] = 0
field_family_profiles_abun[,setdiff(inFungus_cazyme,colnames(field_family_profiles_abun))] = 0
field_family_profiles_abun[field_family_profiles_abun > 0] = 1

termite_family_profile_wide[,setdiff(colnames(field_family_profiles_abun), inTermite_cazyme)] = 0
termite_family_profile_wide[,setdiff(inFungus_cazyme, inTermite_cazyme)] = 0
termite_family_profile_wide[termite_family_profile_wide > 0] = 2

termitomyces_family_profiles_wide[,setdiff(colnames(field_family_profiles_abun), inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[,setdiff(inTermite_cazyme, inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[termitomyces_family_profiles_wide > 0] = 3


comb_data_field = rbind(field_family_profiles_abun, termite_family_profile_wide, termitomyces_family_profiles_wide)

comb_data_reorder = as.data.frame(t(comb_data_field)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present)) %>%
  select(!c("PRE", "POST"))

comb_data_reorder_x = comb_data_reorder[order(-rowSums(comb_data_reorder != 0)),]

cols1 = colorRamp2(c(0.5, 1, 2, 3), c("white", "#03387B", "#8FD485", "#F3D84F"))

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 3/complementarity_cazyme_field.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder, col = cols1, cluster = T, track.height = 0.5)
dev.off()
circos.clear()
#### Field timepoint plot end ####


#### All timepoint plot ####
#==========================#

# Adding ECs not present in the gut microbiome be ARE present in the termite and fungus to the dataframe and setting them to 0
allCazyme_profiles_abun = family_profiles_abun %>%
  dcast( Timepoint ~ family, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "Timepoint")

allCazyme_profiles_abun[,setdiff(inTermite_cazyme,colnames(allCazyme_profiles_abun))] = 0
allCazyme_profiles_abun[,setdiff(inFungus_cazyme,colnames(allCazyme_profiles_abun))] = 0


# Set field number
allCazyme_profiles_abun[1, ] <- ifelse(allCazyme_profiles_abun[1, ] > 0, 3, 0)
# Set post number
allCazyme_profiles_abun[2, ] <- ifelse(allCazyme_profiles_abun[2, ] > 0, 2, 0)
# Set pre number
allCazyme_profiles_abun[3, ] <- ifelse(allCazyme_profiles_abun[3, ] > 0, 1, 0)

# Adding ECs not present in the termite  be ARE present in the gut microbiome and fungus to the dataframe and setting them to 0
termite_family_profile_wide[,setdiff(colnames(allCazyme_profiles_abun), inTermite_cazyme)] = 0
termite_family_profile_wide[,setdiff(inFungus_cazyme, inTermite_cazyme)] = 0
termite_family_profile_wide[termite_family_profile_wide > 0] = 4 # setting those present in the termite to 2

# Adding ECs not present in the fungus  be ARE present in the gut microbiome and termite to the dataframe and setting them to 0
termitomyces_family_profiles_wide[,setdiff(colnames(allCazyme_profiles_abun), inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[,setdiff(inTermite_cazyme, inFungus_cazyme)] = 0
termitomyces_family_profiles_wide[termitomyces_family_profiles_wide > 0] = 5 # setting those present in the fungus to 3


comb_data_all = rbind(allCazyme_profiles_abun, termite_family_profile_wide, termitomyces_family_profiles_wide)

comb_data_reorder_all = as.data.frame(t(comb_data_all)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present))


cols1 = colorRamp2(c(0.5, 1, 2, 3, 4, 5), c("white", "#55C8FA", "#227EEE", "#03387B", "#8FD485", "#F3D84F"))


pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 3/complementarity_cazyme_all.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder_all, col = cols1, cluster = F, track.height = 0.5)
dev.off()
circos.clear()

#### All timepoint plot end ####


microAll = unique(family_profiles_abun$family)

length(microAll)

sum(inFungus_cazyme %in% microAll)/length(microAll)

sum(inTermite_cazyme %in% microAll)/length(microAll)
