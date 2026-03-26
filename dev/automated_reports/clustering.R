#*******************************************************************
# Title: clustering.R
# Description: Clustering based on SNP data.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Aurélie Fischer

rm(list = objects())

library(stringr)
library(data.table)
library(dplyr)

# Parameters
# Specifying path to output dir
output_dir <- getwd()
# Specifying path to phylogeny results
phylogeny_dir <- "../phylogeny"
# Specifying output file of annotation_report.R
today <- Sys.Date()
annotations <- paste0(today, "_Epitrack_annotation_report.Rdata")
rm(today)
# Specifying at which distance groups from clustering must be built (in terms of SNPs)
chosen_distances <- c(5, 15, 50, 100)

cat(sprintf("Loading %s", annotations))
{
  load(annotations)
  rm(annotations)
}

cat(sprintf("Entering phylogeny directory: %s\n", phylogeny_dir))
{
  setwd(phylogeny_dir)
}

# Get all cluster analysis output files
fnames <- dir(recursive = TRUE, full.names = TRUE)
fnames <- fnames[!(fnames %in% c("./command.sh", "./clusters.tsv"))]
clusters <- unique(str_match(fnames, "cluster_[^/]+/\\w+")) # example: c("cluster_1/paeruginosa")
cat(sprintf("  Phylogeny directory contains %i files in %i clusters.\n", length(fnames), length(clusters)))

# Clustering
cat("Scanning MLDIST file from IQtree output or SNP matrix computed from Snippy output.\n")
{
  # Get snp matrix in output files
  # Here if there's gubbins snp matrix take it otherwise keep snippy snp matrix
  # List all files with minicluster/snippy/snp_matrix.tsv
  # change, add : minicluster/gubbins/core_snp_matrix.tsv
  mldist_fnames <- fnames[grepl("minicluster_[0-9]+/snippy/snp_matrix.tsv", fnames)]
  mldist_fnames <- c(mldist_fnames,
                     fnames[grepl("minicluster_[0-9]+/gubbins/core_snp_matrix.tsv", fnames)])
  # Get minicluster name along with cluster name
  mldist_clusters <- str_match(mldist_fnames, "cluster_[^/]+/\\w+/minicluster_\\d+")
  # Check that there is 1 file for 1 cluster in list
  stopifnot(length(mldist_fnames) == length(mldist_clusters))
  # Name list items with cluster id
  names(mldist_fnames) <- mldist_clusters
  cat(sprintf("  Found %i clusters with a MLDIST report.\n", length(unique(mldist_clusters))))
  # Read snp matrix content and associate it to right cluster name
  mldist_reports <- list()
  for(cluster in clusters) {
    mldist_reports[[ cluster ]] <- list()
    for(minicl in unique(mldist_clusters)){
      if(grepl(cluster, minicl, fixed = TRUE)){
        # Stocker dans mldist_reports, la matrice de snps du fichier donné par mldist_fnames[minicl]
        matrix_files <- mldist_fnames[grepl(minicl, mldist_fnames)]
        matrix_files <- ifelse(TRUE %in% grepl("gubbins/core_snp_matrix.tsv", matrix_files),
                               matrix_files[which(grepl("gubbins/core_snp_matrix.tsv", matrix_files))],
                               matrix_files)
        mldist_reports[[ cluster ]][[ minicl ]] <- fread( matrix_files, skip = 1 )
        setnames(mldist_reports[[ cluster ]][[ minicl ]], names(mldist_reports[[ cluster ]][[ minicl ]]), c("NAME", mldist_reports[[ cluster ]][[ minicl ]]$V1))
        mldist_reports[[ cluster ]][[ minicl ]] <- mldist_reports[[ cluster ]][[ minicl ]][,-1]
        if(TRUE %in% grepl("Reference", names(mldist_reports[[ cluster ]][[ minicl ]]))){
          # Removing Reference line in snp matrix because we don't work on it
          total_items <- ncol(mldist_reports[[ cluster ]][[ minicl ]])
          mldist_reports[[ cluster ]][[ minicl ]] <- mldist_reports[[ cluster ]][[ minicl ]][1:(total_items-1),1:(total_items-1)]
          
        } else {print("hey")}
      }
    }
  }
  rm(mldist_fnames, mldist_clusters, cluster, minicl, total_items, matrix_files)
}

cat("Transform genetic distance (substitions/sites) into substitions number.\n")
{
  genetic_distances <- list()
  snps_df <- data.table()
  for(cluster in names(mldist_reports)) {
    main_cl <- mldist_reports[[ cluster ]]
    genetic_distances[[ cluster ]] <- list()
    for(minicluster in names(main_cl)) {
      mini_cl <- main_cl[[ minicluster ]]
      # Transform info into distance 
      genetic_distances[[ cluster ]][[ minicluster ]] <- as.dist(mini_cl)
      # Add number of SNPs if it's cluster with only 2 samples
      if ((length(names(main_cl))==1) & (length(names(mini_cl)) == 2)) {
        mini_snps_df <- data.table(genome = names(mini_cl),
                                   cluster_id = rep(str_extract(cluster, "cluster_[0-9]+")))
        mini_snps_df$nb_snps <- genetic_distances[[ cluster ]][[ minicluster ]]
        snps_df <- rbind.data.frame(snps_df, mini_snps_df)
      }
    }
  }
  rm(mldist_reports, cluster, main_cl, mini_cl, minicluster, mini_snps_df)
}

cat("Loading outliers.\n")
{
  # All outliers are listed in a file named no_clonal_samples.tsv
  outliers_fnames <- fnames[grepl("no_clonal_samples.tsv", fnames)]
  outliers_clusters <- str_match(outliers_fnames, "cluster_[^/]+/\\w+")
  names(outliers_fnames) <- outliers_clusters
  cat(sprintf("  Found %i clusters with outliers detected.\n", length(outliers_clusters)))
  # Set them in a list with cluster id
  outliers_reports <- list()
  for(cluster in clusters) {
    print(cluster)
    outliers_reports[[ cluster ]] <- list()
    outlier_data <- read.table(outliers_fnames[cluster], sep = "\t", header = TRUE)
    outliers_reports[[ cluster ]] <- outlier_data$sample_names
  }
  rm(outlier_data, outliers_clusters, cluster, outliers_fnames)
}

cat("Executing clustering.\n")
{
  clustering <- list()
  clustering_df <- data.frame()
  # Main cluster scale [Example: cluster_2/ecloacae]
  for(cluster in names(genetic_distances)) {
    main_cl <- genetic_distances[[ cluster ]]
    subclustering_df <- data.table()
    print(paste0("Analysing big ", cluster))
    # Minicluster scale [Example: cluster_2/ecloacae - minicluster_1]
    for(minicluster in names(main_cl)) {
      mini_cl <- main_cl[[ minicluster ]]
      num_mini_cl <- str_extract(minicluster, "\\d+$")
      print(paste0("Minicluster n° ", num_mini_cl))
      # Do hierarchical clustering
      cl <- hclust(mini_cl, method = "single")
      # Cut groups for each distance threshold chosen before
      distances_df <- data.table(genome = cl$labels,
                                 cluster_id = rep(str_extract(cluster, "cluster_[0-9]+"), length(cl$labels)))
      # Computing clustering for each chosen threshold
      for(dist in chosen_distances){
        groups <- cutree(cl, h = dist)
        # Store data in a dataframe for this specific distance
        tmp <- data.table(names(groups), str_extract(cluster, "cluster_[0-9]+"), paste0(num_mini_cl, '-', groups))
        colnames(tmp) <- c("genome", "cluster_id", paste0("clust_dist_", dist))
        # Merge results in a dataframe combining all distances of clustering
        distances_df <- merge(distances_df, tmp)
        }
      # Merge results for this minicluster with those from other miniclusters
      print("Is this rbindlist ok ? [Line 158]")
      subclustering_df <- rbind(subclustering_df, distances_df)
      print("Yep")
    }
    # Appending outliers to clustering results
    for(outlier in outliers_reports[[cluster]]){
      # Create random number 
      out_index <- which(outliers_reports[[cluster]] == outlier)
      outlier_line <- data.table(outlier,  
                                 rep(str_extract(cluster, "cluster_[0-9]+")), 
                                 paste0("out-", out_index), paste0("out-", out_index), paste0("out-", out_index), paste0("out-", out_index))
      colnames(outlier_line) <- c("genome", "cluster_id", paste0("clust_dist_", chosen_distances) )
      print("Is this rbindlist OK ? [Line 170]")
      subclustering_df <- rbind(subclustering_df, outlier_line)
      print("Yep")
    }
    # Now let's sort results, to get biggest cluster names 'A' -> smallest cluster 'Z'
    for(dist in chosen_distances){
      print(paste0("Distance analysée = ", dist, " SNPs"))
      # Alphabet in R = LETTERS
      # If there are more than 26 clusters, name with a second letter AA, AB, AC, ...
      # Set biggest group as cluster 1, etc. to get groups in ascending order
      dist_column <- paste0("clust_dist_",dist)
      sorted <- sort(table(subclustering_df[[paste0("clust_dist_",dist)]]), decreasing = TRUE)
      list_of_letters <- c(LETTERS, paste0("A", LETTERS))
      for(i in seq(1,length(sorted))){
        old_name <- names(sorted)[i]
        index_to_change <- which(subclustering_df[[ dist_column ]] == old_name)
        new_name <- list_of_letters[i]
        subclustering_df[[ dist_column ]][index_to_change] <- new_name
      }
    }
    clustering[[ cluster ]] <- subclustering_df
    print("Is this rbindlist OK ? [Line 191]")
    clustering_df <- rbind(clustering_df, subclustering_df)
    print("Yep")
  }
  rm(cl, groups, sorted, i, cluster, genetic_distances, distances_df, main_cl, subclustering_df,
     tmp, dist, dist_column, index_to_change, mini_cl, minicluster, new_name, num_mini_cl, old_name,
     outlier_line, outliers_reports, out_index, outlier)
}

cat("Merging results with metadata table.\n")
{
  cat("Moving to output directory.\n")
  setwd(output_dir)
  if(nrow(snps_df)!=0){
    clustering_df <- merge(clustering_df, snps_df, key = "genome", all = TRUE)
  }
  clustering_df$cluster_id <- str_extract(clustering_df$cluster_id, "[0-9]+")
  # Remove the old cluster id (reports) and keep cluster id used for analysis (clustering_df)
  reports$cluster_id <- c()
  final_data <- merge(clustering_df, reports)
  today_date <- Sys.Date()
  save(final_data, file = paste0(today_date, "_clone_reports.Rdata"))
  
  rm(clustering_df, reports, today_date, fnames, snps_df)
}

cat("Cleaning environment.")
{
  rm(clustering, final_data, chosen_distances, clusters, output_dir, phylogeny_dir, list_of_letters)
}
