#*******************************************************************
#*GENOME CLUSTERING ANALYSIS
#*V0.1 2023-08-16
#*
#* Clustering based on SNP data
#* 
#* Annotate a subcluster for each sample

#TODO
#How do we set the cut to define a clone ? [currently : arbitrary 15 SNPs]

rm(list = objects())

library(stringr)
library(data.table)
library(dplyr)
# Parameters
# Specifying path to output dir
output_dir <- getwd()
# Specifying path to phylogeny results
phylogeny_dir <- "../phylogeny"
# Specifying output file of annotation_report_2.R
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

# Analyze clusters
fnames <- dir(recursive = TRUE, full.names = TRUE)
clusters <- unique(str_match(fnames, "cluster_\\d+/\\w+")) # example: c("cluster_1/paeruginosa")

cat(sprintf("  Phylogeny directory contains %i files in %i clusters.\n", length(fnames), length(clusters)))

# Clustering
cat("Scanning MLDIST file from IQtree output or SNP matrix computed from Snippy output.\n")
{
  mldist_fnames <- fnames[grepl("snp_matrix.tsv", fnames)]
  mldist_fnames <- mldist_fnames[grepl("snippy", mldist_fnames)]
  mldist_clusters <- str_match(mldist_fnames, "cluster_\\d+/\\w+/minicluster_\\d+")
  stopifnot(length(mldist_fnames) == length(mldist_clusters))
  names(mldist_fnames) <- mldist_clusters
  cat(sprintf("  Found %i clusters with a MLDIST report.\n", length(mldist_clusters)))
  
  mldist_reports <- list()
  for(cluster in clusters) {
    mldist_reports[[ cluster ]] <- list()
    for(minicl in mldist_clusters){
      if(grepl(cluster, minicl, fixed = TRUE)){
        mldist_reports[[ cluster ]][[ minicl ]] <- fread( mldist_fnames[minicl], skip = 1 )
        setnames(mldist_reports[[ cluster ]][[ minicl ]], names(mldist_reports[[ cluster ]][[ minicl ]]), c("NAME", mldist_reports[[ cluster ]][[ minicl ]]$V1))
        mldist_reports[[ cluster ]][[ minicl ]] <- mldist_reports[[ cluster ]][[ minicl ]][,-1]
        # Removing Reference line in snp matrix because we don't work on it
        total_items <- ncol(mldist_reports[[ cluster ]][[ minicl ]])
        mldist_reports[[ cluster ]][[ minicl ]] <- mldist_reports[[ cluster ]][[ minicl ]][1:(total_items-1),1:(total_items-1)]
      }
    }
  }
  rm(mldist_fnames, mldist_clusters, cluster, minicl, total_items)
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
      genetic_distances[[ cluster ]][[ minicluster ]] <- as.dist(mini_cl)
      # Add number of SNPs if it's cluster with only 2 samples
      if ((length(names(main_cl))==1) & (length(names(mini_cl)) == 2)) {
        mini_snps_df <- data.table(genome = names(mini_cl))
        mini_snps_df$nb_snps <- genetic_distances[[ cluster ]][[ minicluster ]]
        snps_df <- rbind.data.frame(snps_df, mini_snps_df)
      }
    }
  }
  rm(mldist_reports, cluster, main_cl, mini_cl, minicluster, mini_snps_df)
}

cat("Loading outliers.\n")
{
  outliers_fnames <- fnames[grepl("no_clonal_samples.tsv", fnames)]
  outliers_clusters <- str_match(outliers_fnames, "cluster_\\d+/\\w+")
  names(outliers_fnames) <- outliers_clusters
  cat(sprintf("  Found %i clusters with outliers detected.\n", length(outliers_clusters)))
  outliers_reports <- list()
  for(cluster in clusters) {
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
  # Main cluster scale
  for(cluster in names(genetic_distances)) {
    main_cl <- genetic_distances[[ cluster ]]
    subclustering_df <- data.table()
    print(paste0("Analysing big ", cluster))
    # Minicluster scale
    for(minicluster in names(main_cl)) {
      mini_cl <- main_cl[[ minicluster ]]
      num_mini_cl <- str_extract(minicluster, "\\d$")
      print(paste0("Minicluster n° ", num_mini_cl))
      # Do hierarchical clustering
      cl <- hclust(mini_cl, method = "single")
      # Cut groups for each distance chosen before
      distances_df <- data.table(genome = cl$labels) 
      # Computing clustering for each genetic distance chosen
      for(dist in chosen_distances){
        groups <- cutree(cl, h = dist)
        # Store data in a dataframe for this specific distance
        tmp <- data.table(names(groups), paste0(num_mini_cl, '-', groups))
        colnames(tmp) <- c("genome", paste0("clust_dist_", dist))
        # Merge results in a dataframe combining all distances of clustering
        distances_df <- merge(distances_df, tmp, key = "genome")
      }
      # Merge results for this minicluster with those from other miniclusters
      subclustering_df <- rbind(subclustering_df, distances_df)
    }
    # Appending outliers to clustering results
    for(outlier in outliers_reports[[cluster]]){
      # Get index of outliers
      out_index <- which(outliers_reports[[cluster]] == outlier)
      outlier_line <- data.table(outlier, paste0("out-", out_index), paste0("out-", out_index), paste0("out-", out_index), paste0("out-", out_index))
      colnames(outlier_line) <- c("genome", paste0("clust_dist_", chosen_distances) )
      subclustering_df <- rbind(subclustering_df, outlier_line)
    }
    # Now let's sort results, to get biggest cluster names 'A' -> smallest cluster 'Z'
    for(dist in chosen_distances){
      print(paste0("Distance analysée = ", dist, " SNPs"))
      # Alphabet in R = LETTERS
      # Set biggest group as cluster 1, etc. to get groups in ascending order
      dist_column <- paste0("clust_dist_",dist)
      sorted <- sort(table(subclustering_df[[paste0("clust_dist_",dist)]]), decreasing = TRUE)
      for(i in seq(1,length(sorted))){
        old_name <- names(sorted)[i]
        index_to_change <- which(subclustering_df[[ dist_column ]] == old_name)
        new_name <- LETTERS[i]
        subclustering_df[[ dist_column ]][index_to_change] <- new_name
      }
    }
    clustering[[ cluster ]] <- subclustering_df
    clustering_df <- rbind(clustering_df, subclustering_df)
  }
  rm(cl, groups, sorted, i, cluster, genetic_distances, distances_df, main_cl, subclustering_df,
     tmp, dist, dist_column, index_to_change, mini_cl, minicluster, new_name, num_mini_cl, old_name,
     outlier_line, outliers_reports, out_index, outlier)
}

cat("Merging results with metadata table.\n")
{
  cat("Moving to output directory.\n")
  setwd(output_dir)
  
  clustering_df <- merge(clustering_df, snps_df, key = "genome", all = TRUE)
  final_data <- merge(clustering_df, reports, key = "genome")
  today_date <- Sys.Date()
  save(final_data, file = paste0(today_date, "_clone_reports.Rdata"))
  
  rm(clustering_df, reports, today_date, fnames, snps_df)
}

cat("Cleaning environment.")
{
  rm(clustering, final_data, chosen_distances, clusters, output_dir, phylogeny_dir)
}
