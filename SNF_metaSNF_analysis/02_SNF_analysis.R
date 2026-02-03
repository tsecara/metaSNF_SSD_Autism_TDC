# Check if the argument C is passed; if not, set default value
args <- commandArgs(trailingOnly = TRUE)
C <- ifelse(length(args) > 0, as.integer(args[1]), 3)  # Default to 3 if not provided

print(paste("Number of Clusters (C):", C))

#######Install packages
setwd(".../SNF_output1/")

custom_lib <- "/projects/tsecara/metaSNF/code/R_packages"
if (!dir.exists(custom_lib)) dir.create(custom_lib, recursive = TRUE)
.libPaths(c(custom_lib, .libPaths()))

#ensure_pkgs <- function(pkgs) {
#  to_install <- pkgs[!vapply(pkgs, function(p)
#    nzchar(system.file(package = p, lib.loc = custom_lib)), logical(1))]
#  if (length(to_install)) {
#    install.packages(to_install, lib = custom_lib,
#                     repos = "https://cloud.r-project.org/",
#                     dependencies = TRUE)
#  }
#}

# Core SNFtool deps that often trip people up
#ensure_pkgs(c(
#  "SNFtool",
#  "ExPosition",
#  "prettyGraphs",
#  "alluvial",
#  "igraph",
#  "clue"         # pulls in 'cluster' of a compatible version
#))

#install.packages("remotes")
#remotes::install_version("SNFtool", version = "2.3.1", repos = "https://cloud.r-project.org/")
#install.packages("devtools")
#devtools::install_github("maxconway/SNFtool")

#install.packages(c("ExPosition","prettyGraphs"), lib = custom_lib, repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(ExPosition,lib.loc=custom_lib)
library(prettyGraphs,lib.loc=custom_lib)
library(SNFtool, lib.loc=custom_lib)

#if (!requireNamespace("remotes", quietly = TRUE))
#  install.packages("remotes", lib = custom_lib, repos = "https://cloud.r-project.org/")

#remotes::install_version("mnormt", version = "2.1.1",
#                         lib = custom_lib, repos = "https://cloud.r-project.org/",
#                         upgrade = "never")

#remotes::install_version("psych", version = "2.4.6.26", repos = "https://cloud.r-project.org/")
library(psych,lib.loc=custom_lib)

#remotes::install_version("dunn.test", version = "1.3.6", repos = "https://cloud.r-project.org/")
library(dunn.test,lib.loc=custom_lib)

#remotes::install_version("ggplot2", version = "3.5.1", repos = "https://cloud.r-project.org/")
library(ggplot2,lib.loc=custom_lib)

#remotes::install_version("dplyr", version = "1.1.4", repos = "https://cloud.r-project.org/")
library(dplyr,lib.loc=custom_lib)

#remotes::install_version("tidyverse", version = "1.3.1", repos = "https://cloud.r-project.org/")
library(tidyverse,lib.loc=custom_lib)

#remotes::install_version("broom", version = "1.0.6", repos = "https://cloud.r-project.org/")
library(broom,lib.loc=custom_lib)

#remotes::install_version("fossil", version = "0.4.0", repos = "https://cloud.r-project.org/")
library(fossil,lib.loc=custom_lib)

#remotes::install_version("corrplot", version = "0.95", repos = "https://cloud.r-project.org/")
library(corrplot,lib.loc=custom_lib)

#remotes::install_version("cluster", version = "2.1.6", repos = "https://cloud.r-project.org/")
library(cluster,lib.loc=custom_lib)

#remotes::install_version("MASS", version = "7.3-61", repos = "https://cloud.r-project.org/")
library(MASS,lib.loc=custom_lib)

#remotes::install_version("ggsignif", version = "0.6.4", repos = "https://cloud.r-project.org/")
library(ggsignif,lib.loc=custom_lib)

#install.packages("mclust",lib = custom_lib, repos = "https://cloud.r-project.org/", dependencies = TRUE)
# (optional belt-and-suspenders, in case CRAN metadata misses a dep)
#for (p in c("mvtnorm","Rcpp")) {
#  if (!nzchar(system.file(package = p, lib.loc = custom_lib))) {
#    install.packages(p, lib = custom_lib, repos = "https://cloud.r-project.org/")
#  }
#}

library(mclust, lib.loc = custom_lib)
#install.packages("fpc")
library(fpc, lib.loc = custom_lib)

# Force R to print errors to .err file
options(error = function() { traceback(2); quit(save="no", status=1, runLast=FALSE) })

# source code for function to determine data integration and clustering across re-sampling using SNF and spectral clustering
source("/projects/tsecara/metaSNF/code/SNF_func/robust_core-clustering_function.R")

#Only use this if you want to recode clusters 
recode_clusters <- function(x) {
  map <- c(
    "1" = 1,
    "2" = 3,
    "3" = 2,
    "4" = 5,
    "5" = 4
  )
  x_chr <- as.character(x)
  out <- map[x_chr]
  # if C < 5, this will just return NA for missing labels, so fix that:
  out[is.na(out)] <- x[is.na(out)]
  as.integer(out)
}

print("####################PACKAGES INSTALLED####################")

#Load in data

#Neurocombatted data FINAL Sample 
subjects <- read.csv(".../full_sample.csv", header=TRUE) #list of subject ideas 
CT <- read.csv(".../CT.csv", header=TRUE)
SA <- read.csv(".../SA.csv", header=TRUE)
VOL <- read.csv(".../VOL.csv", header=TRUE)
FA <- read.csv(".../FA.csv", header=TRUE)
MD <- read.csv(".../MD.csv", header=TRUE)
clinical_scog <- read.csv(".../scog_metaSNF_FINAL2.csv", header=TRUE)
EA <- read.csv(".../EA.csv", header=TRUE)
validation_data <- read.csv(".../demo_combined_reduced_FINAL_FINAL_with_ment.csv", header = TRUE)

#Applying feature reduction using mentalizing score
#New ones 
CT_weight <- read.csv(".../CT_cor_filtered2.csv")
SA_weight <- read.csv(".../SA_cor_filtered2.csv")
VOL_weight <- read.csv(".../VOL_cor_filtered2.csv")
FA_weight <- read.csv(".../FA_cor_filtered2.csv")
MD_weight <- read.csv(".../MD_cor_filtered2.csv")
EA_weight <- read.csv(".../EA_cor_filtered2.csv")

#Neurocombatted data FINAL Sample 
CT <- read.csv(".../combat_data3/CT.csv", header=TRUE)
SA <- read.csv(".../combat_data3/SA.csv", header=TRUE)
VOL <- read.csv(".../combat_data3/VOL.csv", header=TRUE)
FA <- read.csv(".../combat_data3/FA.csv", header=TRUE)
MD <- read.csv(".../combat_data3/MD.csv", header=TRUE)
EA <- read.csv(".../combat_data3/EA.csv", header=TRUE)

CT <- CT[, c("participant_id", CT_weight$feature), drop = FALSE]
SA <- SA[, c("participant_id", SA_weight$feature), drop = FALSE]
VOL <- VOL[, c("participant_id", VOL_weight$feature), drop = FALSE]
FA <- FA[, c("participant_id", FA_weight$feature), drop = FALSE]
MD <- MD[, c("participant_id", MD_weight$feature), drop = FALSE]
EA <- EA[, c("participant_id", EA_weight$feature), drop = FALSE]

#Making sure that participants are the same across all datatypes
CT <- CT[(CT$participant_id %in% clinical_scog$participant_id), ]
SA <- SA[(SA$participant_id %in% clinical_scog$participant_id), ]
VOL <- VOL[(VOL$participant_id %in% clinical_scog$participant_id), ]
FA <- FA[(FA$participant_id %in% clinical_scog$participant_id), ]
MD <- MD[(MD$participant_id %in% clinical_scog$participant_id), ]
EA <- EA[(EA$participant_id %in% clinical_scog$participant_id), ]

train_CT_input = CT %>% dplyr::select(-participant_id)
train_SA_input = SA %>%dplyr::select(-participant_id)
train_VOL_input = VOL %>%dplyr::select(-participant_id)
train_FA_input = FA %>%dplyr::select(-participant_id)
train_MD_input = MD %>%dplyr::select(-participant_id)
train_EA_input = EA %>%dplyr::select(-participant_id)
train_clinical_scog_input = clinical_scog %>%dplyr::select(-participant_id)

# normalizing measures within each data type using a function from the SNF package
train_CT_input = standardNormalization(train_CT_input)
train_SA_input = standardNormalization(train_SA_input)
train_VOL_input = standardNormalization(train_VOL_input)
train_FA_input = standardNormalization(train_FA_input)
train_MD_input = standardNormalization(train_MD_input)
train_EA_input = standardNormalization(train_EA_input)
train_clinical_scog_input = standardNormalization(train_clinical_scog_input)

# setting the parameters (finalized after comparisons using 01_metaSNF.R )
K = 48;		# number of neighbors, usually (10~30), Metacluster C
alpha = 0.6;  	# hyperparameter, usually (0.3~0.8) - Metacluster C
t = 10; 	# Number of Iterations, usually (10~20) 

# creating participant distance matrices using euclidean distances
Dist_CT = SNFtool::dist2(as.matrix(train_CT_input),as.matrix(train_CT_input));
Dist_SA = SNFtool::dist2(as.matrix(train_SA_input),as.matrix(train_SA_input));
Dist_VOL = SNFtool::dist2(as.matrix(train_VOL_input),as.matrix(train_VOL_input));
Dist_FA = SNFtool::dist2(as.matrix(train_FA_input),as.matrix(train_FA_input));
Dist_MD = SNFtool::dist2(as.matrix(train_MD_input),as.matrix(train_MD_input));
Dist_EA = SNFtool::dist2(as.matrix(train_EA_input),as.matrix(train_EA_input));
Dist_scog = SNFtool::dist2(as.matrix(train_clinical_scog_input),as.matrix(train_clinical_scog_input));

# creating participant affinity matrices for each data type
AM_CT = affinityMatrix(Dist_CT,K,alpha)
AM_SA = affinityMatrix(Dist_SA,K,alpha)
AM_VOL = affinityMatrix(Dist_VOL,K,alpha)
AM_FA = affinityMatrix(Dist_FA,K, alpha)
AM_MD = affinityMatrix(Dist_MD,K, alpha)
AM_EA = affinityMatrix(Dist_EA,K, alpha)
AM_scog = affinityMatrix(Dist_scog,K, alpha)

C = C  # Using the value of C passed as a parameter

# calling function to integrate data types using SNF and cluster participants using spectral clustering 
#across resampling 80% of participants 1000 times
print("####################SNF RUNNING####################")
robust.W = RobustCoreClusteringMatrix(feature.affinity.mat.list = list(AM_CT, AM_SA, AM_VOL, AM_EA, AM_scog),exp.num.samples = 1000, num.clusts = C, seed = 123)
#Two matrices - Dense Core Cluster Matrix and Sparse Core Cluster Matrix
dense <- robust.W[1]
dense <- matrix(unlist(dense), ncol = 335, byrow = TRUE)

# displaying clusters
# Save the heatmap plot directly while running the script
output_filename <- paste0(C, "_train__clust_Heatmap.png")

# Start the PNG device to save the plot
png(output_filename, width = 1200, height = 800, res = 100)

displayClustersWithHeatmap(dense, spectralClustering(dense, C))
#displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))

# Close the PNG device (this saves the plot)
dev.off()
print("####################CLUSTER MAP####################")

# calculating normalized mutual information (NMI) based off of the original data types and the clustering similarity matrix
SNF_NMIScores <- rankFeaturesByNMI(list(train_CT_input, train_SA_input, train_VOL_input, train_EA_input, train_clinical_scog_input), dense)
#SNF_NMIScores <- rankFeaturesByNMI(list(train_CT_input, train_SA_input, train_EA_input, train_clinical_scog_input), dense)
# separating and organizing scores
CT_scores <- as.data.frame(SNF_NMIScores[[1]][1])
SA_scores <- as.data.frame(SNF_NMIScores[[1]][2])
VOL_scores <- as.data.frame(SNF_NMIScores[[1]][3])
EA_scores <- as.data.frame(SNF_NMIScores[[1]][4])
scog_scores <- as.data.frame(SNF_NMIScores[[1]][5])
names(CT_scores) <- c("NMI")
names(SA_scores) <- c("NMI")
names(VOL_scores) <- c("NMI")
names(EA_scores) <- c("NMI")
names(scog_scores) <- c("NMI")

CT_scores <- t(as.data.frame(CT_scores))
colnames(CT_scores) = colnames(train_CT_input)
SA_scores <- t(as.data.frame(SA_scores))
colnames(SA_scores) = colnames(train_SA_input)
VOL_scores <- t(as.data.frame(VOL_scores))
colnames(VOL_scores) = colnames(train_VOL_input)
EA_scores <- t(as.data.frame(EA_scores))
colnames(EA_scores) = colnames(train_EA_input)
scog_scores <- t(as.data.frame(scog_scores))
colnames(scog_scores) = colnames(train_clinical_scog_input)

all_scores <- cbind(CT_scores, SA_scores)
all_scores <- cbind(all_scores, VOL_scores)
all_scores <- cbind(all_scores, EA_scores)
all_scores <- cbind(all_scores, scog_scores)
all_scores <- t(all_scores)

directory <- ".../SNF_output1/"
#saving csv of NMI scores and clustering similarity matrix
write.csv(all_scores, paste0(".../SNF_output1/full_INPUTS_all_scores_", C, "_clust_k48_0.6_1000perms_transposed.csv"))
print("####################NMI SCORES####################")
write.matrix(dense, file = file.path(directory, paste0("full_clustering_matrix_", C, "_clust_k48_0.6_1000perms")))
print("####################SIMILARITY MATRIX####################")

# Find cluster labels of individuals using the robust clustering similarity matrix
robust.groups.df = RobustCoreClusteringClusters(core.clustering.list = robust.W,num.clusts = C,verbose = T)

# Recode the groups here
robust.groups.df$groups <- recode_clusters(robust.groups.df$groups)

cat("Cluster counts for C =", C, ":\n")
print(table(robust.groups.df$groups))

#train_subjects <- subjects[subjects$"participant_id" %in% train_subs$participant_id, ]
clusters <- cbind(subjects, robust.groups.df)
table(clusters$groups)

#Reorganize and re-integrating demographics file with projected cluters 
clusters$id <- NULL #Removing ID row as this is not necessary 

#Re-loading in datatypes that will be used in the final integrated dataframe 
ori_demo <- read.csv(".../demo_combined_reduced_FINAL_FINAL_with_ment.csv", header = TRUE)
CT_anal <- read.csv(".../CT.csv", header=TRUE)
SA_anal <- read.csv(".../SA.csv", header=TRUE)
VOL_anal <- read.csv(".../VOL.csv", header=TRUE)
FA_anal <- read.csv(".../FA.csv", header=TRUE)
MD_anal <- read.csv(".../MD.csv", header=TRUE)
scog_anal <- read.csv(".../scog_metaSNF_FINAL2.csv", header=TRUE)
EA_anal <- read.csv(".../EA.csv", header=TRUE)

#Combing the two together 
compare_df <- inner_join(clusters, ori_demo, by = "participant_id")
compare_df <- inner_join(compare_df, SA_anal)
compare_df <- inner_join(compare_df, CT_anal)
compare_df <- inner_join(compare_df, VOL_anal)
compare_df <- inner_join(compare_df, FA_anal)
compare_df <- inner_join(compare_df, MD_anal)
compare_df <- inner_join(compare_df, scog_anal)
compare_df <- inner_join(compare_df, EA_anal)


write.csv(compare_df, paste0(".../SNF_output1/full_RESULTS_", C, "_clust_k48_0.6.csv"))
print("####################CLUSTER ASSIGNMENT####################")

## calculating silouette width for each participant and the silhouette plot
dissim <- 1 - dense
dissim <- as.matrix(dissim)
clusters$groups <- as.integer(clusters$groups)
sil <- silhouette(clusters$groups, dmatrix = dissim)
# Convert groups column to numeric
clusters$groups <- as.numeric(as.character(clusters$groups))
# Check if there are any non-numeric values
non_numeric <- sum(is.na(clusters$groups))

# saving the silhouette plot
# Dynamically create a color vector based on the value of C
colors <- switch(as.character(C),
                 "2" = c(adjustcolor("#F8766D", alpha.f = 0.5), 
                         adjustcolor("#A3A500", alpha.f = 0.5)),
                 "3" = c(adjustcolor("#F8766D", alpha.f = 0.5), 
                         adjustcolor("#A3A500", alpha.f = 0.5), 
                         adjustcolor("#00BF7D", alpha.f = 0.5)),
                 "4" = c(adjustcolor("#F8766D", alpha.f = 0.5), 
                         adjustcolor("#A3A500", alpha.f = 0.5), 
                         adjustcolor("#00BF7D", alpha.f = 0.5), 
                         adjustcolor("#00B0F6", alpha.f = 0.5)),
                 "5" = c(adjustcolor("lightgreen", alpha.f = 0.5), 
                         adjustcolor("seagreen", alpha.f = 0.5), 
                         adjustcolor("aquamarine2", alpha.f = 0.5), 
                         adjustcolor("lightblue", alpha.f = 0.5), 
                         adjustcolor("blue", alpha.f = 0.5)),
                 "6" = c(adjustcolor("#F8766D", alpha.f = 0.5), 
                         adjustcolor("#A3A500", alpha.f = 0.5), 
                         adjustcolor("#00BF7D", alpha.f = 0.5), 
                         adjustcolor("#00B0F6", alpha.f = 0.5), 
                         adjustcolor("#E76BF3", alpha.f = 0.5), 
                         adjustcolor("#D89000", alpha.f = 0.5)),
                 "7" = c(adjustcolor("#F8766D", alpha.f = 0.5), 
                         adjustcolor("#A3A500", alpha.f = 0.5), 
                         adjustcolor("#00BF7D", alpha.f = 0.5), 
                         adjustcolor("#00B0F6", alpha.f = 0.5), 
                         adjustcolor("#E76BF3", alpha.f = 0.5), 
                         adjustcolor("#D89000", alpha.f = 0.5),
                         adjustcolor("darkviolet", alpha.f = 0.5)),
                 "8" = c(adjustcolor("#F8766D", alpha.f = 0.5), 
                         adjustcolor("#A3A500", alpha.f = 0.5), 
                         adjustcolor("#00BF7D", alpha.f = 0.5), 
                         adjustcolor("#00B0F6", alpha.f = 0.5), 
                         adjustcolor("#E76BF3", alpha.f = 0.5), 
                         adjustcolor("#D89000", alpha.f = 0.5),
                         adjustcolor("darkviolet", alpha.f = 0.5),
                         adjustcolor("blue", alpha.f = 0.5))
)

# Saving the silhouette plot with dynamic colors
name = paste0(C, "full_clust_Silhouette_plot_0.7_24.png")
png(name, width = 1200, height = 800, res = 100)

plot(sil, col = colors, border = NA, lwd = 1)

dev.off()
print("####################SILHOUETTES####################")


cluster_stats <-fpc::cluster.stats(d=dissim, clusters$groups, alt.clustering = NULL, noisecluster = FALSE, silhouette = TRUE,wgap = TRUE, sepindex = TRUE)
cluster_stats <- as.data.frame(unlist(cluster_stats))
write.csv(cluster_stats, file = paste0(directory, "full_Cluster_fit_stats_", C, ".csv"), row.names = TRUE)
print("####################FIT STATISTICS####################")
