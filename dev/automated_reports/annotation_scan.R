#*******************************************************************
# Title: annotation_scan.R
# Description: Scans annotation files to feed the annotation database
# Authors: Jean-Philippe Rasigade, Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Jean-Philippe Rasigade, Aurélie Fischer

rm(list = objects())

library(data.table)
library(stringr)
library(readr)
options(java.parameters = "-Xmx5000m") # Avoid heap overflow in Java, see https://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
library(openxlsx)

# PARAMETERS
# Where to write the results
output_dir <- getwd()
# Where to find the genomes
genome_dir <- "./genomes"
# Include an excel output ? (may be lengthy)
excel_output <- FALSE
# Do environment cleaning at the end ? 
clean_env <- TRUE

cat(sprintf("Entering genome directory: %s\n", genome_dir)) 
{
  setwd(genome_dir)
}

# Analyse file names
fnames <- dir(recursive = TRUE, full.names = TRUE)
inames <- unique(str_match(fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"])
cat(sprintf("  Genome directory contains %i files in %i isolates.\n", length(fnames), length(inames)))

##########################
## GATHERING ANNOTATION ##
##########################

# MOB-SUITE Contig_reports
cat("Scanning MOB_SUITE reports.\n")
{
  # Get all files concerning contigs from output
  contig_report_fnames <- fnames[grepl("contig_report.txt", fnames)]
  contig_report_inames <- str_match(contig_report_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(contig_report_fnames) == length(contig_report_inames))
  names(contig_report_fnames) <- contig_report_inames
  cat(sprintf("  Found %i isolates with a MOB_SUITE CONTIG report.\n", length(contig_report_inames)))
  # Read file content and place it into a list with sample_id
  contig_reports <- list()
  for(iname in contig_report_inames) {
    contig_reports[[ iname ]] <- fread( contig_report_fnames[iname] )[ , sample_id := iname]
  }
  contig_reports <- rbindlist(contig_reports)[, key := paste(sample_id, contig_id, sep = "_")]
  # Change "sample_id" with "genome" in colnames
  setnames(contig_reports, "sample_id", "genome")
  rm(contig_report_fnames, contig_report_inames, iname)
}
cat("  Finished scanning MOB_SUITE CONTIG reports.\n")

# MOB-SUITE mge_reports
{
  # Get all files concerning mge reports 
  mge_fnames <- fnames[grepl("mge.report.txt", fnames)]
  mge_inames <- str_match(mge_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(mge_fnames) == length(mge_inames))
  names(mge_fnames) <- mge_inames
  cat(sprintf("  Found %i isolates with a MOB_SUITE MGE report.\n", length(mge_inames)))
  # Read file content and place it into a list with sample id
  mge_reports <- list()
  for(iname in mge_inames) {
    mge_reports[[ iname ]] <- fread( mge_fnames[iname] )[ , sample_id := iname]
  }
  mge_reports <- rbindlist(mge_reports)[, key := paste(sample_id, contig_id, sep = "_")]
  # Change "sample_id" to "genome" in colnames
  setnames(mge_reports, "sample_id", "genome")
  rm(mge_fnames, mge_inames, iname)
}
cat("  Finished scanning MOB_SUITE MGE reports.\n")

# MOB-SUITE mobtyper_reports
{
  # Get files containing mobtyper results
  mobtyper_fnames <- fnames[grepl("mobtyper_results.txt", fnames)]
  mobtyper_inames <- str_match(mobtyper_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(mobtyper_fnames) == length(mobtyper_inames))
  names(mobtyper_fnames) <- mobtyper_inames
  cat(sprintf("  Found %i isolates with a MOB_SUITE MOBTYPER report.\n", length(mobtyper_inames)))
  # Read file content and associated it to sample id in a list
  mobtyper_reports <- list()
  for(iname in mobtyper_inames) {
    mobtyper_reports[[ iname ]] <- fread( mobtyper_fnames[iname] )[ , sample_id := iname]
  }
  mobtyper_reports <- rbindlist(mobtyper_reports)[, key := sample_id]
  # Change sample_id to genome in colnames
  setnames(mobtyper_reports, "sample_id", "genome")
  rm(mobtyper_fnames, mobtyper_inames, iname)
}
cat("  Finished scanning MOB_SUITE Plasmid reports.\n")

# MLST reports
cat("Scanning MLST reports.\n")
{
  # Get all files containing MLST results
  mlst_fnames <- fnames[grepl("mlst.tsv", fnames)]
  mlst_inames <- str_match(mlst_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  names(mlst_fnames) <- mlst_inames
  cat(sprintf("  Found %i isolates with a MLST report.\n", length(mlst_inames)))
  # Read MLST report and include it to a list with sample id
  mlst_reports <- list()
  for(iname in mlst_inames) {
    mlst_reports[[ iname ]] <- fread( mlst_fnames[iname], header = FALSE )[ , V1 := iname]
  }
  mlst_reports <- rbindlist(mlst_reports, fill = TRUE)
  # Set proper colnames
  setnames(mlst_reports, names(mlst_reports),
           c("genome", "mlst_species", "sequence_type", 
             {
               if( ncol(mlst_reports) > 3 )
                 paste("allele", 1:(ncol(mlst_reports) - 3), sep = "_")
               else NULL
             }
           ))
  rm(mlst_fnames, mlst_inames, iname)
}
cat("  Finished scanning MLST reports.\n")

# AMR reports
cat("Scanning AMRFINDER+ reports.\n")
{
  # Get all files with AMRfinder+ results
  amrfinder_report_fnames <- fnames[grepl("amrfinder.tsv", fnames)]
  amrfinder_report_inames <- str_match(amrfinder_report_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(amrfinder_report_fnames) == length(amrfinder_report_inames))
  names(amrfinder_report_fnames) <- amrfinder_report_inames
  cat(sprintf("  Found %i isolates with an AMRFINDER+ report.\n", length(amrfinder_report_inames)))
  # Read report content and place it into a list associated with sample id
  amrfinder_reports <- list()
  for(iname in amrfinder_report_inames) {
    amrfinder_reports[[ iname ]] <- fread( amrfinder_report_fnames[iname] )
    amrfinder_reports[[ iname ]] <- cbind(genome = iname, amrfinder_reports[[ iname ]])
  }
  amrfinder_reports <- rbindlist(amrfinder_reports, fill = TRUE)[, key := paste(genome, `Contig id`, sep = "_")]
  rm(amrfinder_report_fnames, amrfinder_report_inames, iname)
}
cat("  Finished scanning AMRFinder+ reports.\n")

# Sourmash reports
cat("Scanning SOURMASH reports.\n")
{
  # Get all files with Sourmash output
  sourmash_report_fnames <- fnames[grepl("sourmash.csv", fnames)]
  sourmash_report_inames <- str_match(sourmash_report_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(sourmash_report_fnames) == length(sourmash_report_inames))
  names(sourmash_report_fnames) <- sourmash_report_inames
  cat(sprintf("  Found %i isolates with a SOURMASH report.\n", length(sourmash_report_inames)))
  # Place content into a list with sample id
  sourmash_reports <- list()
  for(iname in sourmash_report_inames) {
    sourmash_reports[[ iname ]] <- fread( sourmash_report_fnames[iname] )[, ID := iname]
  }
  sourmash_reports <- rbindlist(sourmash_reports)
  # Set proper colnames (unified "genome" for sample id)
  setnames(sourmash_reports, "ID", "genome")
  rm(sourmash_report_fnames, sourmash_report_inames, iname)
}
cat("  Finished scanning SOURMASH reports.\n")
rm(fnames, genome_dir, inames)
cat("Leaving genome directory.\n")

##############
## END SCAN ##
##############

cat(sprintf("Entering output directory: %s\n", output_dir))
setwd(output_dir)

# Merge reports into 1 big data.table
combined_amr_report <- merge(amrfinder_reports, contig_reports, by = c("key", "genome"), all.x = TRUE)
combined_amr_report <- merge(combined_amr_report, sourmash_reports, by = "genome", all.x = TRUE)
combined_amr_report <- merge(combined_amr_report, mlst_reports, by = "genome", all.x = TRUE)

# 2025-07-02
# Adding a merge especially for amr results and combined amr results
cat("Merge equivalent columns in AMR annotation results.\n")
{ 
  # Protein identifier == Protein id
  # Gene symbol == Element symbol
  # Sequence name == Element name
  # Element type == Type
  # Element subtype == Subtype
  
  amrfinder_reports$`Protein identifier` <- ifelse(!is.na(amrfinder_reports$`Protein identifier`), amrfinder_reports$`Protein identifier`, amrfinder_reports$`Protein id`)
  amrfinder_reports$`Gene symbol` <- ifelse(!is.na(amrfinder_reports$`Gene symbol`), amrfinder_reports$`Gene symbol`, amrfinder_reports$`Element symbol`)
  amrfinder_reports$`Sequence name` <- ifelse(!is.na(amrfinder_reports$`Sequence name`), amrfinder_reports$`Sequence name`, amrfinder_reports$`Element name`)
  amrfinder_reports$`Element type` <- ifelse(!is.na(amrfinder_reports$`Element type`), amrfinder_reports$`Element type`, amrfinder_reports$Type)
  amrfinder_reports$`Element subtype` <- ifelse(!is.na(amrfinder_reports$`Element subtype`), amrfinder_reports$`Element subtype`, amrfinder_reports$Subtype)
  
  combined_amr_report$`Protein identifier` <- ifelse(!is.na(combined_amr_report$`Protein identifier`), combined_amr_report$`Protein identifier`, combined_amr_report$`Protein id`)
  combined_amr_report$`Gene symbol` <- ifelse(!is.na(combined_amr_report$`Gene symbol`), combined_amr_report$`Gene symbol`, combined_amr_report$`Element symbol`)
  combined_amr_report$`Sequence name` <- ifelse(!is.na(combined_amr_report$`Sequence name`), combined_amr_report$`Sequence name`, combined_amr_report$`Element name`)
  combined_amr_report$`Element type` <- ifelse(!is.na(combined_amr_report$`Element type`), combined_amr_report$`Element type`, combined_amr_report$Type)
  combined_amr_report$`Element subtype` <- ifelse(!is.na(combined_amr_report$`Element subtype`), combined_amr_report$`Element subtype`, combined_amr_report$Subtype)
  
}

# Just watching some results
# Sourmash detected genus
combined_amr_report[ !duplicated(genome) , .N, genus][ order(-N)]
# SOC species
combined_amr_report[ !duplicated(genome) , .N, species][ order(-N)]
# ARGs
combined_amr_report[ , .N, `Gene symbol`][ order(-N) ]
# Plasmid primary cluster ID
combined_amr_report[ !duplicated(paste(genome, primary_cluster_id)), .N, primary_cluster_id][ order(-N) ]

####################
## SAVING RESULTS ##
####################

cat("Saving R data format.\n")
today_date <- Sys.Date()
save(mge_reports, mobtyper_reports, amrfinder_reports, sourmash_reports, contig_reports, mlst_reports, combined_amr_report, file = paste0(today_date, "_annotation_scan.Rdata"))

# Readable Excel output
# Had a problem with library xlsx
# Error : .jcall("RJavaTools", "Z", "hasField", .jcast(x, "java/lang/Object"),  : 
# RcallMethod: cannot determine object class
# Solution : use of openxlsx library instead

if(excel_output == TRUE) {
  cat("Saving Excel data format. Please be patient...\n")
  data_to_save <- list(combined_amr_report, mge_reports, amrfinder_reports, sourmash_reports, contig_reports, mlst_reports)
  sheet_names <- c("amr", "mge", "amrfinder", "sourmash", "contig", "mlst")
  write.xlsx(data_to_save, file = paste0(today_date, "_annotation_report.xlsx"), sheetName = sheet_names, rowNames = FALSE)
}

cat("Annotation scan finished. Exiting.\n")

if(clean_env){
  rm(amrfinder_reports, combined_amr_report, contig_reports, mobtyper_reports, mge_reports, mlst_reports, sourmash_reports, today_date)
}

rm(clean_env, excel_output, output_dir)

