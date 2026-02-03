# Set the working directory to a specific path
setwd("...../metaSNF")

### Checking iterations of parameters and their effect on cluster number and NMI scores
install.packages("dplyr")
library(dplyr)
install.packages("devtools")
library(devtools)
devtools::install_github("BRANCHlab/metasnf@v0.6.6")
library(metasnf)
install.packages("car")
library(car)
install.packages("future") #for parallel processing
library(future)
install.packages("future.apply")
library(future.apply)


#Loading in datatypes to be used in clustering and for external-validation 
CT <- read.csv(".../CT.csv", header=TRUE)
SA <- read.csv(".../SA.csv", header=TRUE)
VOL <- read.csv(".../VOL.csv", header=TRUE)
FA <- read.csv(".../FA.csv", header=TRUE)
MD <- read.csv(".../MD.csv", header=TRUE)
clinical_scog <- read.csv(".../scog_metaSNF_FINAL2.csv", header=TRUE)
clinical_nc <- read.csv(".../ncog_metaSNF_FINAL2.csv", header=TRUE)
EA <- read.csv(".../EA.csv", header=TRUE)

#### Putting everything into a list will help get quicker summaries: 
study_data <- list(
  CT,
  SA,
  VOL,
  FA,
  MD,
  clinical_scog,
  EA
)

# The number of rows in each dataframe:
lapply(study_data, dim)

# Whether or not each dataframe has missing values:
lapply(study_data,
       function(x) {
         any(is.na(x))
         any(x == 0)
       }
)

### Having option for standardization. This allows you to specify standardization during metaSNF. 
my_distance_metrics <- generate_distance_metrics_list(
  continuous_distances = list(
    "standard_norm_euclidean" = sn_euclidean_distance
  ),
  discrete_distances = list(
    "standard_norm_euclidean" = sn_euclidean_distance
  ),
  ordinal_distances = list(
    "standard_norm_euclidean" = sn_euclidean_distance
  ),
  categorical_distances = list(
    "standard_norm_euclidean" = gower_distance
  ),
  mixed_distances = list(
    "standard_norm_euclidean" = gower_distance
  ),
  keep_defaults = FALSE
)

############### TO BETTER MATCH STANDARD SNF #################
# Custom function: standardNormalization + squared Euclidean - ensures closer match to standard SNF 
dist2_sew <- function(df, weights_row = NULL) {
  df <- SNFtool::standardNormalization(df)
  if (!is.null(weights_row)) {
    weights <- diag(weights_row, nrow = length(weights_row))
    weights <- sqrt(weights)
    df <- as.matrix(df) %*% weights
  }
  distance_matrix <- as.matrix(stats::dist(df, method = "euclidean"))
  distance_matrix <- distance_matrix^2
  return(distance_matrix)
}

# Use this in the distance metric list
my_distance_metrics <- generate_distance_metrics_list(
  continuous_distances = list(
    "standard_norm_sqeuclidean" = dist2_sew
  ),
  discrete_distances = list(
    "standard_norm_sqeuclidean" = dist2_sew
  ),
  ordinal_distances = list(
    "standard_norm_sqeuclidean" = dist2_sew
  ),
  categorical_distances = list(
    "standard_norm_sqeuclidean" = gower_distance
  ),
  mixed_distances = list(
    "standard_norm_sqeuclidean" = gower_distance
  ),
  keep_defaults = FALSE
)

#Creates datalist to be included in the metaSNF
data_list2 <- generate_data_list(
  list(CT, "cortical_thickness", "neuroimaging", "continuous"),
  list(SA, "cortical_surface_area", "neuroimaging", "continuous"),
  list(VOL, "subcortical_volume", "neuroimaging", "continuous"),
  list(EA, "empathic_accuracy", "neuroimaging", "continuous"),
  list(clinical_scog, "scog_variable", "demographics", "continuous"),
  uid = "participant_id"
)


#Summarizing the included data in the list
summarize_dl(data_list2)


#Pre-selecting the cluster size, ensuring clusters size of 3-8 is selected. 
my_spectral_eigen <- function(similarity_matrix) {
  estimated_n <- SNFtool::estimateNumberOfClustersGivenGraph(
    W = similarity_matrix,
    NUMC = 3:8
  )
  number_of_clusters <- estimated_n$`Eigen-gap best`
  solution <- SNFtool::spectralClustering(
    similarity_matrix,
    number_of_clusters
  )
  solution_data <- list("solution" = solution, "nclust" = number_of_clusters)
  return(solution_data)
}
my_spectral_rot <- function(similarity_matrix) {
  estimated_n <- SNFtool::estimateNumberOfClustersGivenGraph(
    W = similarity_matrix,
    NUMC = 3:8
  )
  number_of_clusters <- estimated_n$`Rotation cost best`
  solution <- SNFtool::spectralClustering(
    similarity_matrix,
    number_of_clusters
  )
  solution_data <- list("solution" = solution, "nclust" = number_of_clusters)
  return(solution_data)
}
clust_algs_list <- generate_clust_algs_list(
  "my_first_function" = my_spectral_eigen,
  "my_second_function" = my_spectral_rot,
  disable_base = TRUE # important!
)

###### Now going to run the metaSNF to see the difference
#Generating the settings matrix 

#FINAL SAMPLE SNF
set.seed(42)
settings_matrix2 <- generate_settings_matrix(
  data_list2,
  nrow = 289, #643 is the greatest number of potential combinations possible 
  min_k = 20,
  max_k = 50,
  t_values = 10,
  dropout_dist = "none",
  distance_metrics_list = my_distance_metrics,
  clustering_algorithms = clust_algs_list,
  possible_snf_schemes = 1,
  seed = 42
)

######Running SNF for all the rows in the settings matrix 
#With parallel processing 
solutions_matrix2 <- batch_snf(
  data_list2, 
  settings_matrix2,
  distance_metrics_list = my_distance_metrics,
  clust_algs_list = clust_algs_list,
  processes = 24
  )

#You can pull the clustering results out of each row using: 
cluster_solutions2 <- get_cluster_solutions(solutions_matrix2)
head(cluster_solutions2)

#Calculating ARI across clusting solutions 
solutions_matrix_aris <- calc_aris(solutions_matrix2, processes = 24) 
meta_cluster_order <- get_matrix_order(solutions_matrix_aris)

#Installing package required to select meta-clusters 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ComplexHeatmap", "S4Vectors", "IRanges"))
BiocManager::install("InteractiveComplexHeatmap")
library(InteractiveComplexHeatmap)

#Plotting the ARI heatmap 
adjusted_rand_index_heatmap(
  solutions_matrix_aris,
  order = meta_cluster_order
)
shiny_annotator(ari_hm)

#Denoting meta-cluster boundaries 
split_vec <- c(60, 146, 201, 227, 273, 284) #all input features 

#replotting with clear meta-cluster boudaries 
ari_mc_hm <- adjusted_rand_index_heatmap(
  solutions_matrix_aris,
  order = meta_cluster_order,
  split_vector = split_vec
)
ari_mc_hm

#Loading in the data that will be used to extract target variables and for extternal validation 
validation_data <- read.csv("/projects/tsecara/metaSNF/data/demographics_and_cog/combined/demo_combined_reduced_FINAL_FINAL_with_ment.csv", header = TRUE)
validation_data <- validation_data[(validation_data$participant_id %in% CT$participant_id),]

#Extracting target variable(s) for analysis 
BSFS_data <- validation_data[, c("participant_id", "bsfs_total")]
mentalizing_data <- validation_data[, c("participant_id", "mentalizing")]

target_list <- generate_data_list(
  list(BSFS_data, "functional_outcomes", "behaviour", "numeric"),
  list(mentalizing_data, "mentalizing_scog", "behaviour", "numeric"),
  uid = "participant_id"
)

# Only looking at our out-of-model p-values
extended_solutions_matrix <- extend_solutions(
  solutions_matrix2,
  target_list = target_list,
  processes = 28
)

# Re-running to calculate the p-value for every single input and out-of-model
# feature:
extended_solutions_matrix <- extend_solutions(
  solutions_matrix2,
  data_list = data_list2,
  target_list = target_list,
  processes = 28
)

#Provides an annotated plot to show trands across metaclusters 
annotated_ari_hm <- adjusted_rand_index_heatmap(
  solutions_matrix_aris,
  order = meta_cluster_order,
  split_vector = split_vec,
  data = extended_solutions_matrix,
  top_hm = list(
    "Mentalizing p-value" = "mentalizing_pval",
    #"BSFS Total Score p-value" = "bsfs_total_pval",
    "Overall outcomes p-value" = "mean_pval"
  ),
  bottom_bar = list(
    "Number of Clusters" = "nclust",
    "K Hyperparam" = "k",
    "Alpha" = "alpha"
  ),
  annotation_colours = list(
    #"TASIT3 Lies p-value" = colour_scale(
    #  extended_solutions_matrix$"tasit3_lies_val_pval",
    #  min_colour = "green",
    #  max_colour = "black"
    "Mentalizing p-value" = colour_scale(
      extended_solutions_matrix$"mentalizing_pval",
      min_colour = "green",
      max_colour = "black"
    ),
    #"TASIT3 Sarcasm p-value" = colour_scale(
    #  extended_solutions_matrix$"tasit3_sar_val_pval",
    #  min_colour = "purple",
    #  max_colour = "black"
    #),
    #"BSFS Total Score p-value" = colour_scale(
    #  extended_solutions_matrix$"bsfs_total_pval",
    #  min_colour = "blue",
    #  max_colour = "black"
    #),
    "Overall outcomes p-value" = colour_scale(
      extended_solutions_matrix$"mean_pval",
      min_colour = "lightblue",
      max_colour = "black"
    )
  )
)


annotated_ari_hm

rep_solutions <- get_representative_solutions(
  solutions_matrix_aris,
  split_vector = split_vec,
  order = meta_cluster_order,
  extended_solutions_matrix
)

###Instead now choosing only for some specific regions:
# Example to remove columns containing "mrisdp" and "other_string"
rep_solutions_red <- rep_solutions %>%
  dplyr::select(-contains("X17"), #EA task cortical regions 
         #-contains("MD"), 
         #-contains("surf"), #surface area 
         #-contains("_FA"), #FA values 
         #-contains("thick"), #cortical thickness
         #-contains(".rh"), #right hemisphere EA subcortical regions
         #-contains(".lh"), #left hemisphere EA subcortical regions 
         #-contains("vol") #subcortical volume 
         ) 

mc_manhattan2 <- mc_manhattan_plot(
  rep_solutions_red,
  data_list = data_list2,
  target_list = target_list,
  point_size = 4,
  threshold = 0.05, #can change between 0.01 and 0.05
  text_size = 12,
  domain_colours = c(
    "neuroimaging" = "cadetblue",
    "demographics" = "purple",
    "behaviour" = "firebrick"
  )
)
mc_manhattan2

#Once determining the top clustering solution, extracting it and examining it in terms of extended matrix solution 
solutions_matrix_aris2 <- solutions_matrix_aris[meta_cluster_order, meta_cluster_order]

#Selecting the specific meta-cluster of interest
solutions_matrix_aris3 <- solutions_matrix_aris2[1:60, 1:60] # A
solutions_matrix_aris3 <- solutions_matrix_aris2[61:146, 61:146] # B
solutions_matrix_aris3 <- solutions_matrix_aris2[147:201, 147:201] # C
solutions_matrix_aris3 <- solutions_matrix_aris2[201:227, 201:227] # D
solutions_matrix_aris3 <- solutions_matrix_aris2[228:273, 228:273] # E
solutions_matrix_aris3 <- solutions_matrix_aris2[273:284, 273:284] # F

#Extracting the top 10 best solutions:
# Find the indices of the top 10 highest values in the matrix

# Create a copy of the matrix to avoid modifying the original
matrix_mod <- solutions_matrix_aris3

# Replace the diagonal values with NA
diag(matrix_mod) <- NA

# Compute the average ARI score for each row, ignoring NA values
average_ari_scores <- rowMeans(matrix_mod, na.rm = TRUE)

# Find the indices of the top 10 highest average ARI scores
top_10_avg_indices <- order(average_ari_scores, decreasing = TRUE)[1:10]

# Extract the row labels for the top 10 rows
row_labels_with_top_10_avg_scores <- rownames(matrix_mod)[top_10_avg_indices]
row_labels_with_top_10_avg_scores

# Retrieve the corresponding rows from extended_solutions_matrix
final_solutions <- extended_solutions_matrix %>%
  filter(row_id %in% row_labels_with_top_10_avg_scores)

# Turn average_ari_scores into a data frame with row_id + score
ari_df <- data.frame(
  row_id = as.numeric(names(average_ari_scores)),
  average_ari_score = as.numeric(average_ari_scores)
)

# Join with final_solutions by row_id
final_solutions <- final_solutions %>%
  left_join(ari_df, by = "row_id")

# Reorder so that average_ari_score comes right after 'k'
final_solutions <- final_solutions %>%
  relocate(average_ari_score, nclust, mentalizing_pval, .after = k)


# Display the final solutions
print(final_solutions)
write.csv(final_solutions, "/projects/tsecara/metaSNF/final_solutions_ARI.csv")
