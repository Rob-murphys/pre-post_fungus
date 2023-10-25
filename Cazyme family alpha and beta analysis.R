library(vegan)
library(tibble)
library(dplyr)
library(effectsize)
library(pairwiseAdonis)

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


#### Alpha diversity analysis ####
#================================#

### Shannon ###
shannon_df_cazyme = as.data.frame(vegan::diversity(family_profiles_abun_wide, index = "shannon")) %>%
  rename(shannon_cazyme = `vegan::diversity(family_profiles_abun_wide, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata, by = c("rowname" = "Sample")) %>%
  column_to_rownames()

  ## Anova ##
shannon_cazyme.aov = aov(shannon_cazyme ~ Timepoint * Matriline, data = shannon_df_cazyme)
summary(shannon_cazyme.aov)
TukeyHSD(shannon_cazyme.aov, which = "Timepoint")
cohens_f(shannon_cazyme.aov)

  ## Plot it ##
cols = c("PRE" = "#55C8FA", "POST" = "#227EEE", "FIELD" = "#03387B")
shannon_cazyme_plot = ggplot(shannon_df_cazyme, aes(x = factor(Timepoint, level = c("PRE", "POST", "FIELD")), y = shannon_cazyme))+
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

ggsave(plot = shannon_cazyme_plot, filename = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Supplamentary/shannon_cazyme_plot.pdf", width = 7, height = 7)

### Observed richness ###
## 16S ##
obs_df_cazyme = family_profiles_abun_wide %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>%
  estimateR(.) %>%
  as.data.frame() %>%
  dplyr::slice(.data=., 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample.id") %>%
  left_join(metadata, by = c("sample.id" = "Sample")) %>%
  column_to_rownames(var = "sample.id") %>%
  rename(obs_cazyme = S.obs)

  ## Anova ##
obs_cazyme.aov = aov(obs_cazyme ~ Timepoint * Matriline, data = obs_df_cazyme)
summary(obs_cazyme.aov)
TukeyHSD(obs_cazyme.aov, which = "Timepoint")
cohens_f(obs_cazyme.aov)

  ## Plot it ##
obs_cazyme_plot = ggplot(obs_df_cazyme, aes(x = factor(Timepoint, level = c("PRE", "POST", "FIELD")), y = obs_cazyme))+
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

ggsave(plot = obs_cazyme_plot, filename = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Supplamentary/obs_cazyme_plot.pdf", width = 7, height = 7)


#### Alpha diversity analysis end ####


adonis.res = adonis2(family_profiles_abun_wide ~ Timepoint*Matriline, data = metadata, permutations = 9999, method = "bray", by = "terms")

adonis_pairwise.res = pairwise.adonis2(family_profiles_abun_wide ~ Timepoint * Matriline , data = metadata, method = "bray", nperm = 9999)

#### Lignin substrate comparision ####

substrate_profiles_abun_wide = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/substrate_profiles_abun_wide.rds")

lignin = substrate_profiles_abun_wide %>%
  left_join(metadata, by = c("sample_name" = "Sample")) %>%
  column_to_rownames(var="sample_name") %>%
  select(c("lignin", "Timepoint"))


res.aov = aov(lignin ~ Timepoint, data = lignin)
summary(res.aov)
TukeyHSD(res.aov)
