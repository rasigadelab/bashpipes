# Title: preparation.R
# Description: Prepares a CSV file listing all FASTA file of plasmids carrying an ARG (either blaOXA-48, any blaNDM or any blaVIM) inside a cluster.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Aurélie Fischer

rm(list = ls())

library(data.table)
library(stringr)

# List samples in the analyses based on directories present in "genomes" repo
path_to_results <- "Plasmid_results"
sample_list <- dir(paste0(path_to_results, "/genomes"))

# Search for information around plasmids
reports <- "2025-04-30_Epitrack_annotation_report.Rdata"
scan <- "2025-04-30_annotation_scan.Rdata"
load(reports)
load(scan)
# Filter samples of interest
sample_info <- reports[genome %in% sample_list,c("genome", "args_primary_clusters")]
sample_info_scan <- mobtyper_reports[genome %in% sample_list,]
rm(reports, sample_list, sourmash_reports, mobtyper_reports, amrfinder_reports, combined_amr_report,
   contig_reports, mge_reports, mlst_reports, scan)

# Get blaOXA-48, blaNDM, blaVIM arg
path_to_results_nf <- gsub("~", "Documents", path_to_results)
arg <- c("blaOXA-48", "blaNDM", "blaVIM")
path_to_arg_df <- data.table()
for(s in sample_info$genome){
  tmp <- sample_info[genome==s,]
  for(g in arg){
    # If one of blaOXA-48;blaVIM or blaNDM is found in carried ARGs
    if(grepl(g, tmp$args_primary_clusters)){
      if(!(g == "blaOXA-48")){
        # Get complete gene name for blaVIM or blaNDM
        genes <- unlist(str_extract_all(tmp$args_primary_clusters, paste0(g, "-[0-9]+")))
      } else {
        genes <- g
      }
      for(gene in genes){
        # Get plasmid info (MobSuite ID + Inc + Path to FASTA)
        plasmid_id <- str_extract(tmp$args_primary_clusters, paste0("[^,]+:",gene))
        plasmid_id <- gsub(paste0(":", gene), "", plasmid_id)
        plasmid_inc <- sample_info_scan[genome == s & primary_cluster_id == plasmid_id,c("rep_type(s)")]
        plasmid_inc <- gsub("/", "-", plasmid_inc)
        path_to_plasmid <- paste0(path_to_results_nf, "/genomes/", s, "/mob_recon/plasmid_", plasmid_id, ".fasta")
        # Prepare new row containing (sample, ARG, plasmid carrying arg, inc type, path to FASTA)
        new_row <- cbind(s, gene, plasmid_id, plasmid_inc, path_to_plasmid)
        # Add new row to final df that lists all plasmid-ARG present in samples.
        path_to_arg_df <- rbind(path_to_arg_df, new_row)
      }
    }
  }

}
rm(sample_info, tmp, arg, g, gene, path_to_plasmid, plasmid_id, s, new_row, sample_info_scan,
   plasmid_inc)

# Output to tsv file
fwrite(path_to_arg_df, paste0(path_to_results, "/plasmid_locations.tsv"), sep = "\t", col.names = FALSE)
