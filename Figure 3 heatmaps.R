setwd("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output")
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(tibble)
library(dplyr)

#### Read in data ####
#====================#
## Metadata ##
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/sample_metadata.csv", stringsAsFactors = FALSE)

  ### Cazyme family ###

## Metagenomes ##
family_profiles_abun_wide = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/family_profiles_abun_wide.rds")

## Termitomyces ##
termitomyces_family_profiles = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_family_profiles.RDS")

## M. natalensis ##
termite_family_profile = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_family_profiles.RDS")

  ### Cazyme substrate ###

## Metagenomes ##
substrate_profiles_abun_wide = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/substrate_profiles_abun_wide.rds")

## Termitomyces ##
termitomyces_substrate_profiles = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/termitomyces_substrate_profiles.RDS")

## M. natalensis ##
termite_substrate_profile = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/dbcan4_output/mNatalensis_substrate_profiles.RDS")


#### Read in data end ####

#### Cazyme family heatmap ####
#=============================#

### Preparing datafrmes for heatmap ###
termitomyces_family_profiles$fungus_present = 1

microXfungus = intersect(colnames(family_profiles_abun_wide),termitomyces_family_profiles$family)
termitomyces_fams = data.frame(family = microXfungus, fungus_present = rep("present", length(microXfungus)))

termite_family_profile$termite_present = 2

microXtermite = intersect(colnames(family_profiles_abun_wide),termite_family_profile$family)
termite_fams = data.frame(family = microXtermite, termite_present = rep("present", length(microXtermite)))

family_profiles_abun_wide = family_profiles_abun_wide %>% 
  left_join(metadata, by = c("sample_name" = "Sample")) %>% 
  arrange(desc(Timepoint)) %>%
  dplyr::select(-c(Timepoint, Matriline)) %>%
  column_to_rownames(var = "sample_name")

### Heatmap annotations ###
annotations.cols = as.data.frame(t(family_profiles_abun_wide)) %>% 
  rownames_to_column(var = "rn") %>% 
  select("rn") %>%
  left_join(termitomyces_fams, by = c("rn" = "family")) %>% 
  left_join(termite_fams, by = c("rn" = "family")) %>%
  column_to_rownames(var = "rn") %>%
  replace(is.na(.), "not_present")

annotations.rows = metadata %>% 
  dplyr::select(Sample, Timepoint) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Sample") %>%
  arrange(desc(Timepoint))
annotations.rows$Timepoint = factor(annotations.rows$Timepoint, levels = c("PRE", "POST", "FIELD"))

colours = list(
  Timepoint = c(PRE = "#55C8FA", POST = "#227EEE", FIELD = "#03387B"),
  termite_present = c(present = "#94C785", not_present = "white"),
  fungus_present = c(present = "#F4D94F", not_present = "white"))                               

## log10 transfrm data ##
family_profiles_abun_wide_log10 = log10(family_profiles_abun_wide)

family_profiles_abun_wide_log10[family_profiles_abun_wide_log10 == "-Inf"] = -(1-as.numeric(range(family_profiles_abun_wide_log10[family_profiles_abun_wide_log10 != "-Inf"])[1]))# removing "inf" spawned from log10()

## Plot ##
pdf(file = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 3/cazyme_family_heatmap.pdf", width = 25, height = 10)
ComplexHeatmap::pheatmap(family_profiles_abun_wide_log10, fontsize = 16, 
                         annotation_legend = FALSE,
                         annotation_col = annotations.cols,
                         annotation_colors = colours,
                         show_colnames = FALSE, 
                         annotation_names_row = FALSE,
                         row_split = annotations.rows$Timepoint,
                         cluster_row_slices = FALSE)
dev.off()

#### Cazyme family heatmap end ####



#### Cazyme substrate heatmap ####
#===============================#

termitomyces_substrate_profiles$fungus_present = 1

microXfungus = intersect(colnames(substrate_profiles_abun_wide),termitomyces_substrate_profiles$Substrate)
termitomyces_subs = data.frame(substrate = microXfungus, fungus_present = rep("present", length(microXfungus)))

termite_substrate_profile$termite_present = 2

microXtermite = intersect(colnames(substrate_profiles_abun_wide),termite_substrate_profile$Substrate)
termite_subs = data.frame(substrate = microXtermite, termite_present = rep("present", length(microXtermite)))

substrate_profiles_abun_wide = substrate_profiles_abun_wide %>% 
  left_join(metadata, by = c("sample_name" = "Sample")) %>% 
  arrange(desc(Timepoint)) %>%
  dplyr::select(-c(Timepoint, Matriline)) %>%
  column_to_rownames(var = "sample_name")

### Heatmap annotations ###
annotations.cols = as.data.frame(t(substrate_profiles_abun_wide)) %>% 
  rownames_to_column(var = "rn") %>% 
  select("rn") %>%
  left_join(termitomyces_subs, by = c("rn" = "substrate")) %>% 
  left_join(termite_subs, by = c("rn" = "substrate")) %>%
  column_to_rownames(var = "rn") %>%
  replace(is.na(.), "not_present")

annotations.rows = metadata %>% 
  dplyr::select(Sample, Timepoint) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Sample") %>%
  arrange(desc(Timepoint))
annotations.rows$Timepoint = factor(annotations.rows$Timepoint, levels = c("PRE", "POST", "FIELD"))

colours = list(
  Timepoint = c(PRE = "#55C8FA", POST = "#227EEE", FIELD = "#03387B"),
  termite_present = c(present = "#94C785", not_present = "white"),
  fungus_present = c(present = "#F4D94F", not_present = "white"))                               

## log10 transfrm data ##
substrate_profiles_abun_wide_log10 = log10(substrate_profiles_abun_wide)

substrate_profiles_abun_wide_log10[substrate_profiles_abun_wide_log10 == "-Inf"] = -(1-as.numeric(range(substrate_profiles_abun_wide_log10[substrate_profiles_abun_wide_log10 != "-Inf"])[1]))# removing "inf" spawned from log10()
## column names rotation
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("ComplexHeatmap::pheatmap")
)
## Plot ##
pdf(file = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/Figures/Figure 3/cazyme_substrate_heatmap.pdf", width = 25, height = 10)
ComplexHeatmap::pheatmap(substrate_profiles_abun_wide_log10, fontsize = 16, 
                         annotation_legend = FALSE,
                         annotation_col = annotations.cols,
                         annotation_colors = colours,
                         show_colnames = TRUE,
                         angle_col = "45",
                         annotation_names_row = FALSE,
                         row_split = annotations.rows$Timepoint,
                         cluster_row_slices = FALSE,)
dev.off()
