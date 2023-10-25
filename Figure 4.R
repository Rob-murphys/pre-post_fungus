setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")

library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)

## Loading in the data and formatting it ##
#=========================================#

## Metadata ##
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/sample_metadata.csv", stringsAsFactors = FALSE)
metadata$Matriline = as.factor(metadata$Matriline)
metadata$Timepoint = as.factor(metadata$Timepoint)
metadata = metadata[order(metadata$Sample),] #ordering alphabetically

## Cazyme family data ##
family_profiles_abun_wide = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/family_profiles_abun_wide.rds") %>%
  column_to_rownames(var = "sample_name")

## AZCL data ##

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

#### PCOA of cazyme family data ####
#==================================#

cols = c("PRE" = "#55C8FA", "POST" = "#227EEE", "FIELD" = "#03387B")

cazyme.dist = vegdist(family_profiles_abun_wide, "bray")
cazyme.pcoa = pcoa(cazyme.dist)
cazyme.pcoa.scores = as.data.frame(cazyme.pcoa$vectors[,1:2])

cazyme.pcoa.scores = cazyme.pcoa.scores %>% 
  rownames_to_column(var = "rn")  %>% 
  inner_join(metadata, by = c("rn" = "Sample")) %>%
  column_to_rownames(var = "rn")


PCOA_1_2.cent = aggregate(cbind(Axis.1, Axis.2) ~ Timepoint, data = cazyme.pcoa.scores, FUN = mean)

PCOA_1_2.segs <- merge(cazyme.pcoa.scores, setNames(PCOA_1_2.cent, c('Timepoint','segAxis.1','segAxis.2')),
                                     by = 'Timepoint', sort = FALSE)

cazyme_PCOA_1_2.time.matri = ggplot(cazyme.pcoa.scores, aes(x = Axis.1, y = Axis.2))+
  geom_point(size = 7, aes(colour = Timepoint, shape = Matriline))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  theme(axis.text.y = element_text(colour = "black", size = 32), 
        axis.text.x = element_text(colour = "black", size = 32), 
        legend.position = "none", 
        axis.title.y = element_text(size = 32), 
        axis.title.x = element_text(size = 32, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = PCOA_1_2.segs, mapping = aes(xend = segAxis.1, yend = segAxis.2, colour = Timepoint),
               alpha = 0.2, linewidth = 1)+
  geom_point(data = PCOA_1_2.cent, size = 1)+
  geom_text(data = PCOA_1_2.cent, aes(label = Timepoint), size = 12)+
  xlab(paste("PCo 1 (", round(cazyme.pcoa$values$Relative_eig[1], digits = 3)*100, "%)", sep=""))+
  ylab (paste("PCo 2 (", round(cazyme.pcoa$values$Relative_eig[2], digits = 3)*100, "%)", sep=""))

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 4/pannel A1.pdf", cazyme_PCOA_1_2.time.matri, width = 12, height = 12)

#### PCOA of cazyme family data end ####


#### PCOA of AZCL data ####
#=========================#

cols = c("PRE" = "#55C8FA", "POST" = "#227EEE", "FIELD" = "#03387B")

AZCL.dist = vegdist(clean_azcl_data_wide, "bray")
AZCL.pcoa = pcoa(AZCL.dist)
AZCL.pcoa.scores = as.data.frame(AZCL.pcoa$vectors[,1:2])

AZCL.pcoa.scores = AZCL.pcoa.scores %>% 
  rownames_to_column(var = "rn")  %>% 
  inner_join(azcl_metadata, by = c("rn" = "Sample_ID")) %>%
  column_to_rownames(var = "rn")


AZCL_PCOA_1_2.cent = aggregate(cbind(Axis.1, Axis.2) ~ Timepoint, data = AZCL.pcoa.scores, FUN = mean)

AZCL_PCOA_1_2.segs <- merge(AZCL.pcoa.scores, setNames(AZCL_PCOA_1_2.cent, c('Timepoint','segAxis.1','segAxis.2')),
                       by = 'Timepoint', sort = FALSE)

AZCL_PCOA_1_2.time.matri = ggplot(AZCL.pcoa.scores, aes(x = Axis.1, y = Axis.2))+
  geom_point(size = 7, aes(colour = Timepoint))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  theme(axis.text.y = element_text(colour = "black", size = 32), 
        axis.text.x = element_text(colour = "black", size = 32), 
        legend.position = "none", 
        axis.title.y = element_text(size = 32), 
        axis.title.x = element_text(size = 32, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = AZCL_PCOA_1_2.segs, mapping = aes(xend = segAxis.1, yend = segAxis.2, colour = Timepoint),
               alpha = 0.2, linewidth = 1)+
  geom_point(data = AZCL_PCOA_1_2.cent, size = 1)+
  geom_text(data = AZCL_PCOA_1_2.cent, aes(label = Timepoint), size = 12)+
  xlab(paste("PCo 1 (", round(AZCL.pcoa$values$Relative_eig[1], digits = 3)*100, "%)", sep=""))+
  ylab (paste("PCo 2 (", round(AZCL.pcoa$values$Relative_eig[2], digits = 3)*100, "%)", sep=""))

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 4/pannel A2.pdf", AZCL_PCOA_1_2.time.matri, width = 12, height = 12)

#### PCOA of AZCL data end ####


#### Boxplot of AZCL data ####
#============================#

## Parsing AZCL metadata ##
#subs_rm = c("barley_beta_glucan", "chitosan", "dextran", "casein", "collagen") # substraights to remove
subs_rm = c("barley_beta_glucan","casein", "collagen")

azcl_meta = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/AZCL_data/azcl_enzymes.csv", stringsAsFactors = FALSE)

clean_azcl_meta = azcl_meta %>%filter(!AZCL_Substrate %in% subs_rm)

## Parsing data ##

azcl_data = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/AZCL_data/azcl_data.csv", stringsAsFactors = TRUE)


clean_azcl_data = azcl_data %>% 
  filter(Caste1 != "soldier") %>% 
  filter(!AZCL_Substrate %in% subs_rm) # removing specified substrates and soldier caste



clean_azcl_data_meta = clean_azcl_data %>%
  left_join(clean_azcl_meta[colnames(clean_azcl_meta) %in% c("AZCL_Substrate", "Enzyme")], by = "AZCL_Substrate")


# clean_azcl_data_meta$EC = as.factor(clean_azcl_data_meta$EC)
# clean_azcl_data_meta$EC = ordered(clean_azcl_data_meta$EC, levels = c("3.2.1.1", "3.2.1.99", "3.2.1.78", "3.2.1.4", 
#                                                                       "3.2.1.151", "3.2.1.8", "3.2.1.89", 
#                                                                       "3.2.1.39"))


clean_azcl_data_meta$Enzyme.f = factor(clean_azcl_data_meta$Enzyme.y, levels = c("endo-dextranase", "endo-1,3-beta-glucanase",
                                                                              "endo-cellulase", "endo-chitosanase","endo-1,4-beta-galactanase", 
                                                                              "endo-alpha-1,5-L-arabinanase","endo-beta-1,4-xylanase", 
                                                                              "endo-1,4-beta-xyloglucanase", "endo-1,4-beta-xylanase",
                                                                              "endo-1,4-beta-D-mannanase", "alpha-amylase"))

clean_azcl_data_meta$Timepoint= ordered(clean_azcl_data_meta$Timepoint, levels = c("PRE", "POST", "FIELD")) 

cols = c("PRE" = "#55C8FA", "POST" = "#227EEE", "FIELD" = "#03387B")

azcl_plot = ggplot(clean_azcl_data_meta, aes(x = factor(Timepoint, level = c("PRE", "POST", "FIELD")), 
                                             y = Averages_technical_replicates, colour = Timepoint))+
  geom_boxplot()+
  geom_jitter()+
  scale_x_discrete(labels = c("Pre", "Post", "Field"))+
  scale_colour_manual(values = cols)+
  facet_wrap(~Enzyme.f, scales = "fixed", nrow = 2)+
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "none", 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Timepoint")+
  ylab ("Enzymatic activity")

ggsave(plot = azcl_plot, filename = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 4/pannel Bv2.pdf", width = 35, height = 15, units = "cm")

#### Boxplot of AZCL data end ####