library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

#### Import data ####
#===================#
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/metadata.csv")

metadata_noAncistro = metadata %>%
  filter(Species == "Macrotermes") %>%
  filter(BGI_seq_QC == "OK")

## Unormalised ##
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

#### Compute plotting mectrics ####
#=================================#

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
chao_df_16S$Timepoint = factor(chao_df_16S$Timepoint, levels = c("Pre", "Post", "Field"))
chao_df_16S$Caste = factor(chao_df_16S$Caste, levels = c("worker", "soldier", "fungus"))
#### Compute plotting mectrics end ####

#### Plots ####
cols = c("Pre" = "#55C8FA", "Post" = "#227EEE", "Field" = "#03387B")

## Pannel A ##
ddpcr_plot_16S = ggplot(ddpcr_data_filt, aes(x = Timepoint, y = X16S_per_GutFungusSample, fill = Timepoint, shape = Caste, colour = Timepoint))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5)+
  geom_beeswarm(cex = 2, size = 3)+
  scale_y_log10()+
  #coord_cartesian(ylim = c(1e4, 1.5e8))+
  #scale_y_continuous(breaks = c(1e4, 1e5, 1e6, 1e7, 1e8))
  facet_wrap(~Caste, nrow = 1, scales = "fixed")+
  ylab("16S copies")+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none", axis.title.y = element_text(size = 18), 
        panel.spacing = unit(3, "lines"))
  

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 2/pannel A.pdf", ddpcr_plot_16S, width = 20, height = 8)

## Pannel B
chao1_plot_16S = ggplot(chao_df_16S, aes(x = Timepoint, y = chao1_16S, fill = Timepoint, shape = Caste, colour = Timepoint))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5)+
  geom_beeswarm(cex = 2, size = 3)+
  facet_wrap(~Caste, scales = "fixed", nrow = 1)+
  xlab("Timepoint")+
  ylab("Chao1 richness")+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme_pubr()+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"),
        panel.spacing = unit(3, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 2/pannel B.pdf", chao1_plot_16S, width = 20, height = 8)

combined_plot = ggarrange(ddpcr_plot_16S, chao1_plot_16S, nrow = 2,align = "v")

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 2/combined_plot.pdf", combined_plot, width = 15, height = 8)

