library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(vegan)
library(ape)
library(ggpubr)

#### Read in data ####
#====================#

## Metadata ##
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/sample_metadata.csv", stringsAsFactors = FALSE)
ncyc_metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_pathways.csv")[,-3]# loading in ncyc gene metadata without the 3rd column

### ncyc data ###

## Metagenomes ##
ncyc_profiles_abun_wide = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun_wide.rds") %>%
  column_to_rownames(var = "sample_name")

## Termitomyces ##
termitomyces_ncyc_profiles = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/termitomyces/ncyc_profiles_termito.txt", header = TRUE)

## M. natalensis ##
termite_ncyc_profile = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/mnatalensis/ncyc_profiles_mNata.txt", header = TRUE)

#### Read in data end ####

#### ncyc heatmap ####
#====================#

### Preparing datafrmes for heatmap ###
termitomyces_ncyc_profiles$fungus_present = 1

microXfungus = intersect(colnames(ncyc_profiles_abun_wide),termitomyces_ncyc_profiles$gene_name)
termitomyces_genes = data.frame(gene_name = microXfungus, fungus_present = rep("present", length(microXfungus)))

termite_ncyc_profile$termite_present = 2

microXtermite = intersect(colnames(ncyc_profiles_abun_wide),termite_ncyc_profile$gene_name)
termite_genes = data.frame(gene_name = microXtermite, termite_present = rep("present", length(microXtermite)))

ncyc_profiles_abun_wide_orderd = ncyc_profiles_abun_wide %>%
  rownames_to_column(var = "sample_name") %>%
  left_join(metadata, by = c("sample_name" = "Sample")) %>%
  arrange(desc(Timepoint)) %>%
  dplyr::select(-c(Timepoint, Matriline)) %>%
  column_to_rownames(var = "sample_name")

#ncyc_profiles_abun_wide = column_to_rownames(ncyc_profiles_abun_wide, var = "sample_name")
### Heatmap annotations ###
annotations.cols = as.data.frame(t(ncyc_profiles_abun_wide_orderd)) %>% 
  rownames_to_column(var = "rn") %>% 
  select("rn") %>%
  left_join(termitomyces_genes, by = c("rn" = "gene_name")) %>% 
  left_join(termite_genes, by = c("rn" = "gene_name")) %>%
  left_join(ncyc_metadata, by = c("rn" = "gene_name")) %>%
  column_to_rownames(var = "rn") %>%
  replace(is.na(.), "not_present")

annotations.rows = metadata %>% 
  dplyr::select(Sample, Timepoint) %>%
  arrange(desc(Timepoint), Sample) %>%
  remove_rownames() %>% 
  column_to_rownames(var = "Sample")
annotations.rows$Timepoint = factor(annotations.rows$Timepoint, levels = c("PRE", "POST", "FIELD"))

#ncyc_profiles_abun_wide_orderd = ncyc_profiles_abun_wide[order(match(rownames(ncyc_profiles_abun_wide), rownames(annotations.rows))),]


## Annotation colours ##
colours = list(
  #Timepoint = c(PRE = "#55C8FA", POST = "#227EEE", FIELD = "#03387B"),
  pathway1 = c("Nitrification" = "#B59343", "Denitrification" = "#9C5CA1", 
              "Assimilatory nitrate reduction" = "#75AB3D", "Dissimilatory nitrate reduction" = "#6B7DBD",
              "Nitrogen fixation" = "#CC643E",
              "Anammox" = "#55A87C",
              "Organic degradation and synthesis" = "#CA4569",
              "Others" = "#C575A4"),
  Timepoint = c(PRE = "#55C8FA", POST = "#227EEE", FIELD = "#03387B"),
  termite_present = c(present = "#94C785", not_present = "white"),
  fungus_present = c(present = "#F4D94F", not_present = "white")
)

## log10 transfrm data ##
ncyc_profiles_abun_wide_log10 = log10(ncyc_profiles_abun_wide_orderd)

ncyc_profiles_abun_wide_log10[ncyc_profiles_abun_wide_log10 == "-Inf"] = -(1-as.numeric(range(ncyc_profiles_abun_wide_log10[ncyc_profiles_abun_wide_log10 != "-Inf"])[1]))# removing "inf" spawned from log10()

## Plot ##

pdf(file = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/nitrogen_cycling_heatmap.pdf", width = 25, height = 10)
ComplexHeatmap::pheatmap(ncyc_profiles_abun_wide_log10, fontsize = 16, 
                         annotation_legend = FALSE,
                         annotation_col = annotations.cols,
                         annotation_colors = colours,
                         show_colnames = TRUE, 
                         annotation_names_row = FALSE,
                         row_split = annotations.rows$Timepoint,
                         cluster_row_slices = FALSE, clustering_method = "ward.D2")
dev.off()
#### ncyc heatmap end ####

#### PCOA of ncyc data ####
#=========================#
cols = c("PRE" = "#55C8FA", "POST" = "#227EEE", "FIELD" = "#03387B")

ncyc.dist = vegdist(ncyc_profiles_abun_wide, "bray")
ncyc.pcoa = pcoa(ncyc.dist)
ncyc.pcoa.scores = as.data.frame(ncyc.pcoa$vectors[,1:2])

ncyc.pcoa.scores = ncyc.pcoa.scores %>% 
  rownames_to_column(var = "rn")  %>% 
  inner_join(metadata, by = c("rn" = "Sample")) %>%
  column_to_rownames(var = "rn")


PCOA_1_2.cent = aggregate(cbind(Axis.1, Axis.2) ~ Timepoint, data = ncyc.pcoa.scores, FUN = mean)

PCOA_1_2.segs <- merge(ncyc.pcoa.scores, setNames(PCOA_1_2.cent, c('Timepoint','segAxis.1','segAxis.2')),
                       by = 'Timepoint', sort = FALSE)

ncyc_PCOA_1_2.time.matri = ggplot(ncyc.pcoa.scores, aes(x = Axis.1, y = Axis.2))+
  geom_point(size = 7, aes(colour = Timepoint, shape = Matriline))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  theme(axis.text.y = element_text(colour = "black", size = 32), 
        axis.text.x = element_text(colour = "black", size = 32), 
        legend.position = "none", axis.title.y = element_text(size = 32), 
        axis.title.x = element_text(size = 32, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = PCOA_1_2.segs, mapping = aes(xend = segAxis.1, yend = segAxis.2, colour = Timepoint),
               alpha = 0.2, linewidth = 1)+
  geom_point(data = PCOA_1_2.cent, size = 1)+
  geom_text(data = PCOA_1_2.cent, aes(label = Timepoint), size = 12)+
  xlab(paste("PCo 1 (", round(ncyc.pcoa$values$Relative_eig[1], digits = 3)*100, "%)", sep=""))+
  ylab (paste("PCo 2 (", round(ncyc.pcoa$values$Relative_eig[2], digits = 3)*100, "%)", sep=""))

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/pannel B.pdf", ncyc_PCOA_1_2.time.matri, width = 12, height = 12)

#### PCOA of ncyc data end ####


#### ncyc complementary ####
#==========================#

## Read in data ##

## Metadata ##
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/sample_metadata.csv", stringsAsFactors = FALSE)
ncyc_metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_pathways.csv")[,-3]# loading in ncyc gene metadata without the 3rd column

### ncyc data ###

## Metagenomes ##
ncyc_profiles_abun = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun.rds")  %>%
  mutate(present = ifelse(rel_abun > 0, 1, 0)) %>%
  left_join(metadata, by = c("sample_name" = "Sample")) %>%
  select(-c("gene_id", "Contig", "mean_cov")) %>%
  left_join(ncyc_metadata, by = c("gene_name" = "gene_name"))

## Termitomyces ##
termitomyces_ncyc_profiles = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/termitomyces/ncyc_profiles_termito.txt", header = TRUE)

termitomyces_ncyc_profiles_wide = as.data.frame(t(termitomyces_ncyc_profiles %>%
                                                    column_to_rownames(var = "gene_name") %>%
                                                    mutate(present = 1))) %>%
                                                    slice(5)
inFungus_ncyc = colnames(termitomyces_ncyc_profiles_wide)

## M. natalensis ##
termite_ncyc_profile = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/mnatalensis/ncyc_profiles_mNata.txt", header = TRUE)

termite_ncyc_profile_wide = as.data.frame(t(termite_ncyc_profile %>%
                                              column_to_rownames(var = "gene_name") %>%
                                              mutate(mNatalensis = 1))) %>%
                                              filter(!row.names(.) == "Mnat_gene_v1.2.pep")
inTermite_ncyc = colnames(termite_ncyc_profile_wide)

#======================#
## Pre timepoint plot ##
#======================#

# Adding ECs not present in the gut microbiome be ARE present in the termite and fungus to the dataframe and setting them to 0
preNcyc_profiles_abun = ncyc_profiles_abun %>%
  dcast( Timepoint ~ gene_name, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "Timepoint")

preNcyc_profiles_abun[,setdiff(inTermite_ncyc,colnames(preNcyc_profiles_abun))] = 0
preNcyc_profiles_abun[,setdiff(inFungus_ncyc,colnames(preNcyc_profiles_abun))] = 0
preNcyc_profiles_abun[preNcyc_profiles_abun > 0] = 1 # setting those present in the micoribome to 1


# Adding ECs not present in the termite  be ARE present in the gut microbiome and fungus to the dataframe and setting them to 0
termite_ncyc_profile_wide[,setdiff(colnames(preNcyc_profiles_abun), inTermite_ncyc)] = 0
termite_ncyc_profile_wide[,setdiff(inFungus_ncyc, inTermite_ncyc)] = 0
termite_ncyc_profile_wide[termite_ncyc_profile_wide > 0] = 4 # setting those present in the termite to 2

# Adding ECs not present in the fungus  be ARE present in the gut microbiome and termite to the dataframe and setting them to 0
termitomyces_ncyc_profiles_wide[,setdiff(colnames(preNcyc_profiles_abun), inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[,setdiff(inTermite_ncyc, inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[termitomyces_ncyc_profiles_wide > 0] = 5 # setting those present in the fungus to 3


comb_data_pre = rbind(preNcyc_profiles_abun, termite_ncyc_profile_wide, termitomyces_ncyc_profiles_wide)

comb_data_reorder_pre = as.data.frame(t(comb_data_pre)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present)) %>%
  select(!c("POST", "FIELD"))


cols1 = colorRamp2(c(0.5, 1, 2, 3, 4, 5), c("white", "#55C8FA", "#227EEE", "#03387B", "#8FD485", "white"))

splits = ncyc_metadata %>%
  filter(gene_name %in% rownames(comb_data_reorder_pre)) 
splits$pathway1 = as.factor(splits$pathway1)
splits_orderd = splits[order(match(splits$gene_name, rownames(comb_data_reorder_pre))),]

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/complementarity_ncyc_pre_label.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder_pre, col = cols1, cluster = T, track.height = 0.5, 
               split = splits_orderd$pathway1, show.sector.labels = T)
dev.off()
circos.clear()

#=======================#
## Post timepoint plot ##
#=======================#

postNcyc_profiles_abun = ncyc_profiles_abun %>%
  dcast( Timepoint ~ gene_name, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "Timepoint")

postNcyc_profiles_abun[,setdiff(inTermite_ncyc,colnames(postNcyc_profiles_abun))] = 0
postNcyc_profiles_abun[,setdiff(inFungus_ncyc,colnames(postNcyc_profiles_abun))] = 0
postNcyc_profiles_abun[postNcyc_profiles_abun > 0] = 1 # setting those postsent in the micoribome to 1


# Adding ECs not postsent in the termite  be ARE postsent in the gut microbiome and fungus to the dataframe and setting them to 0
termite_ncyc_profile_wide[,setdiff(colnames(postNcyc_profiles_abun), inTermite_ncyc)] = 0
termite_ncyc_profile_wide[,setdiff(inFungus_ncyc, inTermite_ncyc)] = 0
termite_ncyc_profile_wide[termite_ncyc_profile_wide > 0] = 2 # setting those postsent in the termite to 2

# Adding ECs not postsent in the fungus  be ARE postsent in the gut microbiome and termite to the dataframe and setting them to 0
termitomyces_ncyc_profiles_wide[,setdiff(colnames(postNcyc_profiles_abun), inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[,setdiff(inTermite_ncyc, inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[termitomyces_ncyc_profiles_wide > 0] = 3 # setting those postsent in the fungus to 3


comb_data_post = rbind(postNcyc_profiles_abun, termite_ncyc_profile_wide, termitomyces_ncyc_profiles_wide)

comb_data_reorder_post = as.data.frame(t(comb_data_post)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present)) %>%
  select(!c("PRE", "FIELD"))

cols1 = colorRamp2(c(0.5, 1, 2, 3), c("white", "#227EEE", "#8FD485", "#F3D84F"))

splits = ncyc_metadata %>%
  filter(gene_name %in% rownames(comb_data_reorder_post)) 
splits$pathway1 = as.factor(splits$pathway1)
splits_orderd = splits[order(match(splits$gene_name, rownames(comb_data_reorder_post))),]

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/complementarity_ncyc_post_label.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder_post, col = cols1, cluster = T, track.height = 0.5, 
               split = splits_orderd$pathway1, show.sector.labels = T)
dev.off()
circos.clear()

#========================#
## Field timepoint plot ##
#========================#

fieldNcyc_profiles_abun = ncyc_profiles_abun %>%
  dcast( Timepoint ~ gene_name, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "Timepoint")

fieldNcyc_profiles_abun[,setdiff(inTermite_ncyc,colnames(fieldNcyc_profiles_abun))] = 0
fieldNcyc_profiles_abun[,setdiff(inFungus_ncyc,colnames(fieldNcyc_profiles_abun))] = 0
fieldNcyc_profiles_abun[fieldNcyc_profiles_abun > 0] = 1 # setting those fieldsent in the micoribome to 1


# Adding ECs not fieldsent in the termite  be ARE fieldsent in the gut microbiome and fungus to the dataframe and setting them to 0
termite_ncyc_profile_wide[,setdiff(colnames(fieldNcyc_profiles_abun), inTermite_ncyc)] = 0
termite_ncyc_profile_wide[,setdiff(inFungus_ncyc, inTermite_ncyc)] = 0
termite_ncyc_profile_wide[termite_ncyc_profile_wide > 0] = 2 # setting those fieldsent in the termite to 2

# Adding ECs not fieldsent in the fungus  be ARE fieldsent in the gut microbiome and termite to the dataframe and setting them to 0
termitomyces_ncyc_profiles_wide[,setdiff(colnames(fieldNcyc_profiles_abun), inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[,setdiff(inTermite_ncyc, inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[termitomyces_ncyc_profiles_wide > 0] = 3 # setting those fieldsent in the fungus to 3


comb_data_field = rbind(fieldNcyc_profiles_abun, termite_ncyc_profile_wide, termitomyces_ncyc_profiles_wide)

comb_data_reorder_field = as.data.frame(t(comb_data_field)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present)) %>%
  select(!c("PRE", "POST"))

cols1 = colorRamp2(c(0.5, 1, 2, 3), c("white", "#03387B", "#8FD485", "#F3D84F"))

splits = ncyc_metadata %>%
  filter(gene_name %in% rownames(comb_data_reorder_field)) 
splits$pathway1 = as.factor(splits$pathway1)
splits_orderd = splits[order(match(splits$gene_name, rownames(comb_data_reorder_field))),]

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/complementarity_ncyc_field_label.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder_field, col = cols1, cluster = T, track.height = 0.5, 
               split = splits_orderd$pathway1, show.sector.labels = T)
dev.off()
circos.clear()

#======================#
## All timepoint plot ##
#======================#

# Adding ECs not present in the gut microbiome be ARE present in the termite and fungus to the dataframe and setting them to 0
allNcyc_profiles_abun = ncyc_profiles_abun %>%
  dcast( Timepoint ~ gene_name, value.var = "present", fun.aggregate = sum) %>%
  column_to_rownames(var = "Timepoint")

allNcyc_profiles_abun[,setdiff(inTermite_ncyc,colnames(allNcyc_profiles_abun))] = 0
allNcyc_profiles_abun[,setdiff(inFungus_ncyc,colnames(allNcyc_profiles_abun))] = 0


# Set field number
allNcyc_profiles_abun[1, ] <- ifelse(allNcyc_profiles_abun[1, ] > 0, 3, 0)
# Set post number
allNcyc_profiles_abun[2, ] <- ifelse(allNcyc_profiles_abun[2, ] > 0, 2, 0)
# Set pre number
allNcyc_profiles_abun[3, ] <- ifelse(allNcyc_profiles_abun[3, ] > 0, 1, 0)

# Adding ECs not present in the termite  be ARE present in the gut microbiome and fungus to the dataframe and setting them to 0
termite_ncyc_profile_wide[,setdiff(colnames(allNcyc_profiles_abun), inTermite_ncyc)] = 0
termite_ncyc_profile_wide[,setdiff(inFungus_ncyc, inTermite_ncyc)] = 0
termite_ncyc_profile_wide[termite_ncyc_profile_wide > 0] = 4 # setting those present in the termite to 2

# Adding ECs not present in the fungus  be ARE present in the gut microbiome and termite to the dataframe and setting them to 0
termitomyces_ncyc_profiles_wide[,setdiff(colnames(allNcyc_profiles_abun), inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[,setdiff(inTermite_ncyc, inFungus_ncyc)] = 0
termitomyces_ncyc_profiles_wide[termitomyces_ncyc_profiles_wide > 0] = 5 # setting those present in the fungus to 3


comb_data_pre = rbind(allNcyc_profiles_abun, termite_ncyc_profile_wide, termitomyces_ncyc_profiles_wide)

comb_data_reorder_all = as.data.frame(t(comb_data_pre)) %>%
  arrange(desc(PRE), desc(POST), desc(FIELD), desc(mNatalensis), desc(present))


cols1 = colorRamp2(c(0.5, 1, 2, 3, 4, 5), c("white", "#55C8FA", "#227EEE", "#03387B", "#8FD485", "#F3D84F"))

splits = ncyc_metadata %>%
  filter(gene_name %in% rownames(comb_data_reorder_all)) 
splits$pathway1 = as.factor(splits$pathway1)
splits_orderd = splits[order(match(splits$gene_name, rownames(comb_data_reorder_all))),]

pdf("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/complementarity_ncyc_all.pdf", height = 9.44, width = 9.44)
circos.heatmap(comb_data_reorder_all, col = cols1, cluster = T, track.height = 0.5, 
               split = splits_orderd$pathway1)
dev.off()
circos.clear()

microAll = unique(ncyc_profiles_abun$gene_name)

length(microAll)

sum(inFungus_ncyc %in% microAll)/length(microAll)

sum(inTermite_ncyc %in% microAll)/length(microAll)

#### ncyc complementary end ####