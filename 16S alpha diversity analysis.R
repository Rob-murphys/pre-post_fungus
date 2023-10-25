library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(effectsize)

#### Import data ####
#===================#
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/metadata.csv")

metadata_noAncistro = metadata %>%
  filter(Species == "Macrotermes") %>%
  filter(BGI_seq_QC == "OK")

pre_post_feature_table = as.data.frame(t(readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/pre_post_feature_table_decontam.rds")))

pre_post_feature_table_filt = pre_post_feature_table %>%
  arrange(as.numeric(row.names(.))) %>%
  filter(row.names(.) %in% metadata_noAncistro$Sample_ID)


ddpcr_data = read.table("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/ddPCR_results_and_metadata.txt", 
                        sep = "\t",
                        header = T)

ddpcr_data_filt = ddpcr_data %>%
  arrange(as.numeric(row.names(.))) %>%
  filter(Species == "Macrotermes") %>%
  filter(Timepoint %in% c("Pre", "Post", "Field"))
ddpcr_data_filt$Timepoint = factor(ddpcr_data_filt$Timepoint, levels = c("Pre", "Post", "Field"))
ddpcr_data_filt$Caste = factor(ddpcr_data_filt$Caste, levels = c("worker", "soldier", "fungus"))
#### Import data end ####

#### Alpha diversity ####
#=======================#

## Shannon ##

shannon_df_16S = as.data.frame(vegan::diversity(pre_post_feature_table_filt, index = "shannon")) %>%
  rename(shannon_16S = `vegan::diversity(pre_post_feature_table_filt, index = "shannon")`) %>%
  rownames_to_column() %>%
  mutate_at(c("rowname"), as.numeric) %>%
  left_join(metadata_noAncistro, by = c("rowname" = "Sample_ID")) %>%
  column_to_rownames()
shannon_df_16S$Timepoint = factor(shannon_df_16S$Timepoint, levels = c("Pre", "Post", "Field"))

  # Anova #
shannon_16S.aov = aov(shannon_16S ~ Timepoint, data = shannon_df_16S)
summary(shannon_16S.aov)
shannon_HSD_16S = TukeyHSD(shannon_16S.aov, which = "Timepoint")
cohens_f(shannon_16S.aov)

## Chao1 ##
chao_df_16S = pre_post_feature_table_filt %>%
  round() %>%
  estimateR(.) %>%
  as.data.frame() %>%
  dplyr::slice(.data=., 2) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate_at(c("rowname"), as.numeric) %>%
  left_join(metadata_noAncistro, by = c("rowname" = "Sample_ID")) %>%
  column_to_rownames() %>%
  rename(chao1_16S = S.chao1)

  # Anova #
chao_16S.aov = aov(chao1_16S ~ Timepoint * Caste , data = chao_df_16S)
summary(chao_16S.aov)
chao_HSD_16S = TukeyHSD(chao_16S.aov)
cohens_f(chao_16S.aov)

#### ddPCR ####
#=============#
ddpcr_data_filt$log_16S_per_sample = log(ddpcr_data_filt$X16S_per_GutFungusSample)

ddpcr_aov.res = aov(ddpcr_data_filt$log_16S_per_sample ~ ddpcr_data_filt$Timepoint * ddpcr_data_filt$Caste)
summary(ddpcr_aov.res)
TukeyHSD(ddpcr_aov.res)
cohens_f(ddpcr_aov.res)

  ### Worker

worker.res = aov(ddpcr_data_filt[ddpcr_data_filt$Caste == "worker",]$log_16S_per_sample ~ ddpcr_data_filt[ddpcr_data_filt$Caste == "worker",]$Timepoint)
summary(worker.res)
TukeyHSD(worker.res)

  ### Soldier
soldier.res = aov(ddpcr_data_filt[ddpcr_data_filt$Caste == "soldier",]$log_16S_per_sample ~ ddpcr_data_filt[ddpcr_data_filt$Caste == "soldier",]$Timepoint)
summary(soldier.res)
TukeyHSD(soldier.res)

  ### Fungus comb
fungus.res = aov(ddpcr_data_filt[ddpcr_data_filt$Caste == "fungus",]$log_16S_per_sample ~ ddpcr_data_filt[ddpcr_data_filt$Caste == "fungus",]$Timepoint)
summary(fungus.res)
TukeyHSD(fungus.res)


#### ddPCR end ####


