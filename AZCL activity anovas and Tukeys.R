library(dplyr)
library(effectsize)
library(vegan)
library(pairwiseAdonis)

#### Analysis of 9 enzymes in Fig 4c ####
#===========================================#

#Set working directory
setwd("/Users/csp839/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Pre-post manuscript revisions/Data analysis AZCL")

#Load metadata
#prior to loading, 'soliders' were manually removed from the dataframe
#prior to loaded the AZCL substrates 'casein', 'collagen', and 'barley-beta-glucan' were manually removed from the dataframe
azcl_data = read.csv("/Users/csp839/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Pre-post manuscript revisions/Data analysis AZCL/AZCL_data.csv", stringsAsFactors = TRUE)
#
#anovas and posthoc test for enzyme activity
res.aov1 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "amylose",])
summary(res.aov1)
cohens_f(res.aov1)
TukeyHSD(res.aov1)

res.aov2 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "arabinan",])
summary(res.aov2)
cohens_f(res.aov2)
TukeyHSD(res.aov2)

res.aov3 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "arabinoxylan",])
summary(res.aov3)
cohens_f(res.aov3)
TukeyHSD(res.aov3)

res.aov4 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "curdlan",])
summary(res.aov4)
cohens_f(res.aov4)
TukeyHSD(res.aov4)

res.aov5 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "galactan",])
summary(res.aov5)
cohens_f(res.aov5)
TukeyHSD(res.aov5)

res.aov6 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "galactan",])
summary(res.aov6)
cohens_f(res.aov6)
TukeyHSD(res.aov6)

res.aov7 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "galactomannan",])
summary(res.aov7)
cohens_f(res.aov7)
TukeyHSD(res.aov7)

res.aov8 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "HE_cellulose",])
summary(res.aov8)
cohens_f(res.aov8)
TukeyHSD(res.aov8)

res.aov9 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "xylan",])
summary(res.aov9)
cohens_f(res.aov9)
TukeyHSD(res.aov9)

res.aov10 = aov(Averages_technical_replicates ~ Timepoint, data = azcl_data[azcl_data$AZCL_Substrate == "xyloglucan",])
summary(res.aov10)
cohens_f(res.aov10)
TukeyHSD(res.aov10)
