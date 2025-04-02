rm(list = ls())

library(data.table)
library(stringr)

# List samples in the analyses
path_to_results <- "~/Projets/35.Plasmina/Plasmid_results"
sample_list <- dir(paste0(path_to_results, "/genomes"))

# Search for information around plasmids
reports <- "prepa/2025-04-02_Epitrack_annotation_report.Rdata"
scan <- "prepa/2025-04-02_annotation_scan.Rdata"
load(reports)
load(scan)

# Filter samples
sample_info <- reports[genome %in% sample_list,c("genome", "args_primary_clusters")]
sample_info_scan <- mobtyper_reports[genome %in% sample_list,]
rm(reports, sample_list, sourmash_reports, mobtyper_reports, amrfinder_reports, combined_amr_report,
   contig_reports, mge_reports, mlst_reports, scan)

# Get blaOXA-48, blaNDM, blaVIM arg
path_to_results_nf <- gsub("~", "/mnt/c/Users/admin/Documents", path_to_results)
arg <- c("blaOXA-48", "blaNDM", "blaVIM")
path_to_arg_df <- data.table()
for(s in sample_info$genome){
  tmp <- sample_info[genome==s,]
  for(g in arg){
    if(grepl(g, tmp$args_primary_clusters)){
      if(!(g == "blaOXA-48")){
        gene <- str_extract(tmp$args_primary_clusters, paste0(g, "-[0-9]+"))
      } else {
        gene <- g
      }
      plasmid_id <- str_extract(tmp$args_primary_clusters, paste0("[^,]+:",gene))
      plasmid_id <- gsub(paste0(":", gene), "", plasmid_id)
      plasmid_inc <- sample_info_scan[genome == s & primary_cluster_id == plasmid_id,c("rep_type(s)")]
      plasmid_inc <- gsub("/", "-", plasmid_inc)
      path_to_plasmid <- paste0(path_to_results_nf, "/genomes/", s, "/mob_recon/plasmid_", plasmid_id, ".fasta")
      # Prepare new row
      new_row <- cbind(s, gene, plasmid_id, plasmid_inc, path_to_plasmid)
      # Add new row to final df
      path_to_arg_df <- rbind(path_to_arg_df, new_row)
    }
  }

}
rm(sample_info, tmp, arg, g, gene, path_to_plasmid, plasmid_id, s, new_row, sample_info_scan,
   plasmid_inc)

# Output to tsv file
fwrite(path_to_arg_df, paste0(path_to_results, "/plasmid_locations.tsv"), sep = "\t", col.names = FALSE)
