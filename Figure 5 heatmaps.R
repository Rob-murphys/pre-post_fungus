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
  rownames_to_column(var = "sample_name")
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
                         cluster_row_slices = FALSE, clustering_method = "ward.D",main = "ward.D")
dev.off()

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
  geom_point(size = 4, aes(colour = Timepoint, shape = Matriline))+
  theme_pubclean()+
  theme_classic()+
  scale_colour_manual(values = cols)+
  theme(axis.text.y = element_text(colour = "black", size = 18), 
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.position = "none", axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  geom_segment(data = PCOA_1_2.segs, mapping = aes(xend = segAxis.1, yend = segAxis.2, colour = Timepoint),
               alpha = 0.2, linewidth = 1)+
  geom_point(data = PCOA_1_2.cent, size = 1)+
  geom_text(data = PCOA_1_2.cent, aes(label = Timepoint), size = 12)+
  xlab(paste("PCo 1 (", round(ncyc.pcoa$values$Relative_eig[1], digits = 3)*100, "%)", sep=""))+
  ylab (paste("PCo 2 (", round(ncyc.pcoa$values$Relative_eig[2], digits = 3)*100, "%)", sep=""))

ggsave("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 5/pannel A.pdf", ncyc_PCOA_1_2.time.matri, width = 12, height = 12)

#### PCOA of ncyc data end ####