rm(list = ls())

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
libraries <- c("argparse", "stringr", "data.table", "dplyr")
for (lib in libraries){usePackage(lib)}

# Add input parameters to script
parser <- ArgumentParser(description='Bacterial Phylogeny - MASH Distance Matrix Analysis', add_help=TRUE)
parser$add_argument("-d", dest="path_to_mash_dist", required=TRUE, help="Path to the folder containing mash distance matrix", nargs='+')
parser$add_argument("-th", dest="threshold", required=TRUE, help="Clustering threshold, in terms of mash distance, to divide cluster into subgroups")
args <- parser$parse_args()

cat("Moving to directory containing mash matrix.\n")
setwd(args$path_to_mash_dist)
getwd()
replicon <- unlist(str_split(getwd(), "/"))
replicon <- replicon[length(replicon)-1]
cat("Loading mash matrix.\n")
{
  # Load file
  x <- scan(file = "mash_dist.tsv", what = 'character', na.strings='NULL', skip = 1, quiet = TRUE)
  # Creating empty square matrix
  dims <- floor(sqrt(length(x) * 2))
  m <- matrix(NA, dims, dims)
  # Loading data in square matrix
  m[upper.tri(m, diag = TRUE)] <- x
  m <- t(m)
  m <- cbind(m, rep(NA, dims))
  # Removing first columns (sample names)
  sample_names <- str_extract(m[,1],"Epi-\\d+")
  m <- m[,-1]
  # Making matrix symetrical
  my_mat_sym <- apply(m, 2, as.numeric)                                        # Duplicate matrix
  my_mat_sym[upper.tri(my_mat_sym)] <- t(my_mat_sym)[upper.tri(my_mat_sym)] # Insert lower to upper matrix
  rm(m, dims, libraries, x)
}
cat("Clustering based on mash matrix.\n")
{
  # Hierarchical clustering based on mash matrix, single linkage method.
  mash_distances <- my_mat_sym
  hc <- hclust(as.dist(mash_distances), method="single")
  mini_clusters <- cutree(hc, h=as.numeric(args$threshold))
  rm(hc, my_mat_sym)
}
cat("Get reference sample for each minicluster.")
{
  mini_cl_df <- data.frame(sample_names, mini_clusters) %>%
    add_count(mini_clusters, name = "cluster_count")
  mash_df <- data.frame(mash_distances)
  colnames(mash_df) <- sample_names
  mini_cl_df$reference <- rep(NA, length(sample_names))
  for(cl in unique(mini_cl_df$mini_clusters)){
    cl_samples <- mini_cl_df[mini_cl_df$mini_clusters==cl,]$sample_names
    if(length(cl_samples) > 1){
      # Take sub-matrix of the cluster
      sub_mash_df <- mash_df[which(sample_names %in% cl_samples), cl_samples]
      # Compute mean column
      mean_cols <- colMeans(sub_mash_df, na.rm=T)
      ref <- names(mean_cols[which.min(mean_cols)])
      mini_cl_df[mini_cl_df$sample_names %in% cl_samples,]$reference <- ref
    }
  }
  # Take reference in biggest cluster
  ref <- mini_cl_df[which.max(mini_cl_df$cluster_count),]$reference
  
  rm(mash_df, mash_distances, sub_mash_df, cl, cl_samples, mean_cols)
}
cat("Output to TSV file.\n")
{
  # If sample alone in a cluster, write it to no_clonal_samples.tsv
  fwrite(mini_cl_df[mini_cl_df$cluster_count == 1,], file = "no_clonal_samples.tsv", sep = "\t", col.names = TRUE)
  # Else, write it to mini_clusters.tsv
  fwrite(mini_cl_df[mini_cl_df$cluster_count > 1,], file = "mini_clusters.tsv", sep = "\t", col.names = TRUE)
  # Write reference in reference_for_phylogeny.tsv
  fwrite(data.table(replicon = replicon, reference = as.character(ref)), file = "reference_for_phylogeny.tsv", sep = "\t", col.names = TRUE)
  rm(mini_cl_df, sample_names, mini_clusters, lib, usePackage, ref)
}
