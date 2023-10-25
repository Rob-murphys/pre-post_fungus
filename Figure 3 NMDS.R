library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ape)
library(ggpubr)

#### Import data ####
#===================#
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/metadata.csv")

metadata_noAncistro = metadata %>%
  filter(Species == "Macrotermes") %>%
  filter(BGI_seq_QC == "OK")

## Normalised ##
pre_post_feature_table.norm = as.data.frame(readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/pre_post_feature_table_normalised.rds"))

pre_post_feature_table_filt.norm = pre_post_feature_table.norm %>%
  arrange(as.numeric(row.names(.))) %>%
  filter(row.names(.) %in% metadata_noAncistro$Sample_ID)

#### Compute plotting mectrics ####
#=================================#

## Bray Curtis ##
bray_dis = vegdist(pre_post_feature_table_filt.norm, method = "bray")
bray_dismat = as.matrix(bray_dis) # as a matrix

#### Compute plotting mectrics end ####


#### NMDS plot ####

cols = c("Pre" = "#55C8FA", "Post" = "#227EEE", "Field" = "#03387B")

pcoa_16S = pcoa(bray_dis)
pcoa_16S.scores = as.data.frame(pcoa_16S$vectors[,1:2])

nmds_16S = vegan::metaMDS(bray_dismat, k = 4, maxit = 999, trymax = 500) # NMDS of the normalized feature table in bray-curtis space
nmds_16S.scores = as.data.frame(nmds_16S$points)

nmds_16S.scores.meta = nmds_16S.scores %>%
  rownames_to_column(var = "rn") %>%
  mutate_at(c("rn"), as.numeric) %>%
  left_join(metadata_noAncistro, by = c("rn" = "Sample_ID")) %>%
  column_to_rownames(var = "rn")

nmds_16S.scores.meta$Timepoint = factor(nmds_16S.scores.meta$Timepoint, levels = c("Pre", "Post", "Field"))



NMDS_1_2.Timepoint.cent = aggregate(cbind(MDS1, MDS2) ~ Timepoint, data = nmds_16S.scores.meta, FUN = mean)

NMDS_1_2.Timepoint.segs <- merge(nmds_16S.scores.meta, setNames(NMDS_1_2.Timepoint.cent, c('Timepoint','oNMDS1','oNMDS2')),
                                 by = 'Timepoint', sort = FALSE)

NMDS_1_2.Timepoint = ggplot(nmds_16S.scores.meta, aes(x = MDS1, y = MDS2))+
  geom_point(size = 4, aes(colour = Timepoint))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  theme(axis.text.y = element_text(colour = "black", size = 32), 
        axis.text.x = element_text(colour = "black", size = 32), 
        legend.position = "none", axis.title.y = element_text(size = 32), 
        axis.title.x = element_text(size = 32, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = NMDS_1_2.Timepoint.segs, mapping = aes(xend = oNMDS1, yend = oNMDS2),
               alpha = 0.2)+
  geom_point(data = NMDS_1_2.Timepoint.cent, size = 1)+
  geom_text(data = NMDS_1_2.Timepoint.cent, aes(label = Timepoint), size = 12)

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 3/16S_NMDS.pdf", NMDS_1_2.Timepoint, width = 12, height = 12)

#### PCOA plot ####
pcoa_16S = pcoa(bray_dis)
pcoa_16S.scores = as.data.frame(pcoa_16S$vectors[,1:2])

pcoa_16S.scores.meta = pcoa_16S.scores %>%
  rownames_to_column(var = "rn") %>%
  mutate_at(c("rn"), as.numeric) %>%
  left_join(metadata_noAncistro, by = c("rn" = "Sample_ID")) %>%
  column_to_rownames(var = "rn")
PCOA_1_2.cent = aggregate(cbind(Axis.1, Axis.2) ~ Timepoint, data = pcoa_16S.scores.meta, FUN = mean)

PCOA_1_2.segs <- merge(pcoa_16S.scores.meta, setNames(PCOA_1_2.cent, c('Timepoint','segAxis.1','segAxis.2')),
                       by = 'Timepoint', sort = FALSE)

PCOA_1_2.Timepoint = ggplot(pcoa_16S.scores.meta, aes(x = Axis.1, y = Axis.2))+
  geom_point(size = 7, aes(colour = Timepoint))+
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
  xlab(paste("PCo 1 (", round(pcoa_16S$values$Relative_eig[1], digits = 3)*100, "%)", sep=""))+
  ylab (paste("PCo 2 (", round(pcoa_16S$values$Relative_eig[2], digits = 3)*100, "%)", sep=""))
ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/figure 3/16S_PCOA.pdf", PCOA_1_2.Timepoint, width = 12, height = 12)

