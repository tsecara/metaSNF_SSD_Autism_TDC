#!/usr/bin/env Rscript

### Script to run resampling stability analysis
library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)
#install.packages("reshape")
library(reshape)

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

subjects <- read.csv(".../full_sample.csv", header=TRUE)
subjects_id <- subjects$participant_id

# Neurocombat-adjusted data
CT  <- read.csv("/projects/tsecara/metaSNF/data/combat_data/combat_data3/CT.csv")
SA  <- read.csv("/projects/tsecara/metaSNF/data/combat_data/combat_data3/SA.csv")
VOL <- read.csv("/projects/tsecara/metaSNF/data/combat_data/combat_data3/VOL.csv")
EA  <- read.csv("/projects/tsecara/metaSNF/data/combat_data/combat_data3/EA.csv")
clinical_scog <- read.csv("/projects/tsecara/metaSNF/data/demographics_and_cog/combined/sep_social_neurocog/final/scog_metaSNF_FINAL2.csv")

# Feature selection weights
CT_weight  <- read.csv("/projects/tsecara/metaSNF/data/weights/correlation_threshold/CT_cor_filtered2.csv")
SA_weight  <- read.csv("/projects/tsecara/metaSNF/data/weights/correlation_threshold/SA_cor_filtered2.csv")
VOL_weight <- read.csv("/projects/tsecara/metaSNF/data/weights/correlation_threshold/VOL_cor_filtered2.csv")
EA_weight  <- read.csv("/projects/tsecara/metaSNF/data/weights/correlation_threshold/EA_cor_filtered2.csv")

# Apply feature reduction
CT  <- CT[,  c("participant_id", CT_weight$feature), drop=FALSE]
SA  <- SA[,  c("participant_id", SA_weight$feature), drop=FALSE]
VOL <- VOL[, c("participant_id", VOL_weight$feature), drop=FALSE]
EA  <- EA[,  c("participant_id", EA_weight$feature), drop=FALSE]

# Keep only participants with cog data
CT  <- CT[CT$participant_id %in% subjects_id, ]
SA  <- SA[SA$participant_id %in% subjects_id, ]
VOL <- VOL[VOL$participant_id %in% subjects_id, ]
EA  <- EA[EA$participant_id %in% subjects_id, ]
clinical_scog  <- clinical_scog[clinical_scog$participant_id %in% subjects_id, ]

#Removing the participant ID
CT  <- CT[ , !(names(CT) %in% "participant_id")]
SA  <- SA[ , !(names(SA) %in% "participant_id")]
VOL <- VOL[ , !(names(VOL) %in% "participant_id")]
EA  <- EA[ , !(names(EA) %in% "participant_id")]
clinical_scog <- clinical_scog[ , !(names(clinical_scog) %in% "participant_id")]

directory <- (".../SNF_output1/stability_testing/")

source(".../bootstrapping2.r")
numboot=1000
nsub=335
K=48
alpha=0.6
t=10
bootsize=0.8
clusters=5

# Setting up which participants will be included in each permutation
permutation_matrix <- bootstrapping_SNF(numboot=numboot, nsub=nsub, bootsize=bootsize)
# Getting clustering solutions for all the permuatations of sampled participants using SNF 
clus_sil <- clustering(perms=permutation_matrix, bootsize=bootsize, K=K, t=t, alpha=alpha, clusters=clusters, CT=CT, SA=SA, VOL=VOL, EA=EA, clinical_scog=clinical_scog)
# Dividing output matrix into clusters and silhouette widths
clus_out <- clus_sil[1:(numboot), ]
silhouette_width <- clus_sil[(numboot+1):(numboot*2), ]

# getting the adjusted rand index between all clustering solutions
list_randindex <- stability(clus_out=clus_out, perms=permutation_matrix) # returns brandindex:adjusted rand index
# Calculate how often each participant is clusted together and the probability that they will be clustered together
percent_agree <- percent_agree(clus_out=clus_out)

write.csv(percent_agree, file=file.path(directory, paste("Percent_agree_5c_1000perms.csv", sep="")))
write.csv(list_randindex, file=file.path(directory, paste("Rand_indices_5c_1000perms.csv", sep="")))

write.csv(clus_out, file=file.path(directory, paste("Adj_rand_indices_5c_1000perms.csv", sep="")))
write.csv(permutation_matrix, file=file.path(directory, paste("permutation_matrix_5c_1000perms.csv", sep="")))
