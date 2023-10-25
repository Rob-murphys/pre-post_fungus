library(dada2)

pre_post_feature_table = readRDS("D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/seqtab.nochim.rds")

silva_taxa = assignTaxonomy(pre_post_feature_table, "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/silva_nr99_v138.1_train_set.fa.gz", 
               multithread=FALSE, verbose=TRUE)

saveRDS(silva_taxa, file = "D:/OneDrive - University of Copenhagen/PhD/Projects/Termite_metagenomes_anvio/16S_work/silva_taxa.RDS")
