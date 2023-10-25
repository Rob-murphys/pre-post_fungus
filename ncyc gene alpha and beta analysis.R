library(tibble)
library(dplyr)
library(vegan)
library(effectsize)
library(ggbeeswarm)
library(ggpubr)
library(ggplot2)
library(pairwiseAdonis)
## Loading in the data and formatting it ##
#=========================================#

## Metadata ##
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/sample_metadata.csv", stringsAsFactors = FALSE)
metadata$Matriline = as.factor(metadata$Matriline)
metadata$Timepoint = as.factor(metadata$Timepoint)
metadata = metadata[order(metadata$Sample),] #ordering alphabetically

## Cazyme family data ##
ncyc_profiles_abun_wide = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/nitrogen_cycling/ncyc_profile_abun_wide.rds") %>%
  column_to_rownames(var = "sample_name")


#### Alpha diversity analysis ####
#================================#

### Shannon ###
shannon_df_ncyc = as.data.frame(vegan::diversity(ncyc_profiles_abun_wide, index = "shannon")) %>%
  rename(shannon_ncyc = `vegan::diversity(ncyc_profiles_abun_wide, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata, by = c("rowname" = "Sample")) %>%
  column_to_rownames()

## Anova ##
shannon_ncyc.aov = aov(shannon_ncyc ~ Timepoint * Matriline, data = shannon_df_ncyc)
summary(shannon_ncyc.aov)
TukeyHSD(shannon_ncyc.aov, which = "Timepoint")
cohens_f(shannon_ncyc.aov)

## Plot it ##
cols = c("PRE" = "#55C8FA", "POST" = "#227EEE", "FIELD" = "#03387B")
shannon_ncyc_plot = ggplot(shannon_df_ncyc, aes(x = factor(Timepoint, level = c("PRE", "POST", "FIELD")), y = shannon_ncyc))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = Timepoint, shape = Matriline), size = 3)+
  xlab("Timepoint")+
  ylab("Shannon diversity")+
  scale_colour_manual(values = cols)+
  scale_x_discrete(labels = c("Pre", "Post", "Field"))+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave(plot = shannon_ncyc_plot, filename = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Supplamentary/shannon_ncyc_plot.pdf", width = 7, height = 7)

### observed richness ###
## 16S ##
obs_df_ncyc = ncyc_profiles_abun_wide %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>%
  estimateR(.) %>%
  as.data.frame() %>%
  dplyr::slice(.data=., 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample.id") %>%
  left_join(metadata, by = c("sample.id" = "Sample")) %>%
  column_to_rownames(var = "sample.id") %>%
  rename(obs_ncyc = S.obs)

## Anova ##
obs_ncyc.aov = aov(obs_ncyc ~ Timepoint * Matriline, data = obs_df_ncyc)
summary(obs_ncyc.aov)
TukeyHSD(obs_ncyc.aov, which = "Timepoint")
cohens_f(obs_ncyc.aov)

## Plot it ##
obs_ncyc_plot = ggplot(obs_df_ncyc, aes(x = factor(Timepoint, level = c("PRE", "POST", "FIELD")), y = obs_ncyc))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 2, aes(colour = Timepoint, shape = Matriline), size = 3)+
  xlab("Timepoint")+
  ylab("Observed richness")+
  scale_colour_manual(values = cols)+
  scale_x_discrete(labels = c("Pre", "Post", "Field"))+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"))

ggsave(plot = obs_ncyc_plot, filename = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Supplamentary/obs_ncyc_plot.pdf", width = 7, height = 7)


#### Alpha diversity analysis end ####

#### Beta diversity analysis ####
adonis.res = adonis2(ncyc_profiles_abun_wide ~ Timepoint*Matriline, data = metadata, permutations = 9999, method = "bray", by = "terms")

adonis_pairwise.res = pairwise.adonis2(ncyc_profiles_abun_wide ~ Timepoint * Matriline , data = metadata, method = "bray", nperm = 9999)
#### Beta diversity analysis end ####