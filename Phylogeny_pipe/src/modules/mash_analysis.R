# Title: mash_analysis.R
# Description: Clustering of genomes using MASH distance.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Aurélie Fischer

rm(list = ls())

# Import/install libraries as needed
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
libraries <- c("optparse", "stringr", "data.table", "dplyr")
for (lib in libraries){usePackage(lib)}

# Add input parameters to script
# So the script could be called from command line
option_list <- list(
  make_option(c("-d", "--path_to_mash_dist"), type="character", help="Path to the folder containing mash distance matrix"),
  make_option(c("-t", "--threshold"), type="double", help="Clustering threshold, in terms of mash distance, to divide cluster into subgroups")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

cat("Moving to directory containing mash matrix.\n")
{
  setwd(args$path_to_mash_dist)
  getwd()
  # Architecture of results directory: /results/phylogeny/replicon/mash
  replicon <- unlist(str_split(getwd(), "/"))
  # Get parent directory to get replicon name
  replicon <- replicon[length(replicon)-1]
}

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
  my_mat_sym <- apply(m, 2, as.numeric) # Duplicate matrix
  my_mat_sym[upper.tri(my_mat_sym)] <- t(my_mat_sym)[upper.tri(my_mat_sym)] # Insert lower to upper matrix
  rm(m, dims, libraries, x)
}

cat("Clustering based on mash matrix.\n")
{
  # Hierarchical clustering based on mash matrix, single linkage method.
  mash_distances <- my_mat_sym
  hc <- hclust(as.dist(mash_distances), method="single")
  # Cut clusters with MASH threshold given in input parameters
  mini_clusters <- cutree(hc, h=as.numeric(args$threshold))
  rm(hc, my_mat_sym)
}

cat("Get reference sample for each minicluster.")
{
  # Create dataframe with samples associated to their minicluster
  mini_cl_df <- data.frame(sample_names, mini_clusters) %>%
    add_count(mini_clusters, name = "cluster_count")
  # Get Mash distances in a df
  mash_df <- data.frame(mash_distances)
  colnames(mash_df) <- sample_names
  # At first no candidate for reference sample
  mini_cl_df$reference <- rep(NA, length(sample_names))
  # Need 1 reference sample for each minicluster
  for(cl in unique(mini_cl_df$mini_clusters)){
    cl_samples <- mini_cl_df[mini_cl_df$mini_clusters==cl,]$sample_names
    if(length(cl_samples) > 1){
      # Take sub-matrix of MASH distances between samples of the minicluster
      sub_mash_df <- mash_df[which(sample_names %in% cl_samples), cl_samples]
      # Compute mean column
      mean_cols <- colMeans(sub_mash_df, na.rm=T)
      # Reference sample is the sample with lowest mean distance with the other samples of minicluster
      ref <- names(mean_cols[which.min(mean_cols)])
      # Attribute the reference to all samples of minicluster
      mini_cl_df[mini_cl_df$sample_names %in% cl_samples,]$reference <- ref
    }
  }
  # Need of 1 reference for phylogeny with all samples from all clusters
  # Take the reference of the biggest minicluster
  ref <- mini_cl_df[which.max(mini_cl_df$cluster_count),]$reference
  if(is.na(ref)){
    # Compute mean column
    mean_cols <- colMeans(mash_df, na.rm=T)
    # If there is no reference as for now, take the sample with lowest mean mash distance with the other samples
    ref <- names(mean_cols[which.min(mean_cols)])
  }
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
