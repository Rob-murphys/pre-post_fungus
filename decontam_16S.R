library(phyloseq)
library(decontam)
library(dplyr)
library(tibble)
source("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/scripts/final scripts/decontam_functions.R")

#### Import data ####
#===================#
metadata = read.csv("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/metadata.csv")

metadata_filt = metadata %>% 
  filter(!Species == "Ancistrotermes") %>%
  filter(!Timepoint %in% c("mock", "unknown label")) %>%
  filter(!BGI_seq_QC == "Weird") %>%
  mutate(Sample_or_Control = case_when(Timepoint == "neg control" ~ "Control Sample", # Create new column with control or sample
                                       TRUE ~ "True Sample"))

pre_post_feature_table = as.data.frame(readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/seqtab.nochim.rds")) 

pre_post_feature_table_filt = pre_post_feature_table %>%
  arrange(as.numeric(row.names(.))) %>%
  filter(row.names(.) %in% metadata_filt$Sample_ID)

taxonomy_table = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/silva_taxa.RDS")
####  Make phyloseq object ####
#=============================#
metadata_phylo = metadata_filt %>%
  column_to_rownames(var = "Sample_ID")

ASV_16S = otu_table(t(pre_post_feature_table_filt), taxa_are_rows = TRUE)
samples_16S = sample_data(metadata_phylo)
taxonomy = tax_table(taxonomy_table)

physeq = phyloseq(ASV_16S, taxonomy, samples_16S)



lib_size(physeq = physeq)

physeq.removed = remove_unwanted(physeq)

# Removing singletons and keeping only those with at least 100 reads in any one sample
physeq.removed.low = remove_low(physeq.removed) # removed 21157 taxa

# Identifying and removing contaminats based on prevalence method
physeq.noncontam = id_contam(physeq.removed.low = physeq.removed.low) # identified 4 contaminants

pre_post_feature_table_clean = otu_table(physeq.noncontam)

saveRDS(pre_post_feature_table_clean, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/pre_post_feature_table_decontam.rds")
