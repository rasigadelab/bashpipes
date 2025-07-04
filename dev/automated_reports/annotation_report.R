#*******************************************************************
#*GENOME ANNOTATION REPORT
#*V0.2 2023-08-02
#*
#* Analyse and report isolate sequencing
#* 
#* Organize tables etc before proceeding

rm(list = objects())

library(data.table)
library(readxl)
library(stringr)
library(lubridate)
library(dplyr)

# Helper function inspired by "curve()" to apply any expression to
# genomes in the 'reports' list

#' @title List apply by object
#' @param lst a list of objects with identical structure
#' @param expr an arbitrary expression to be applied to each element
#' @export
oapply <- function(lst, expr) {
  sexpr <- substitute(expr)
  
  sapply(lst, function(x) eval(sexpr, envir = x, enclos = parent.frame()))
}

# Parameters
# Specify output from annotation_scan.R
today <- Sys.Date()
sample_data <- paste0(today, "_annotation_scan.Rdata")
rm(today)
# Specify metadata file
metadata_file <- "~/Projets/5.Suivi/Metadata_eq_rasigade.xlsx"
# Specify project name : "Epitrack" or "Resistrack"
project <- "Epitrack"
# Specify output directory path
output_dir <- getwd()

cat(sprintf("Loading %s \n", sample_data))
{
  load(sample_data)
  rm(sample_data)
}

cat("Loading metadata. \n")
{
  if(project=="Epitrack"){
    metadata <- data.table(read_excel(metadata_file, sheet = 1, col_types = "text"))
  } else if(project=="Resistrack"){
    metadata <- data.table(read_excel(metadata_file, sheet = 2, col_types = "text"))
  }
  # Set a proper genome key for a better match with sample infos
  setnames(metadata, "SAMPLE_ID", "genome")
  
  rm(metadata_file)
}

cat("Rebuild tables with Epitrack format.\n")
{ 
  raw_table_names <- objects()
  raw <- list()
  for(tbname in raw_table_names){
    item <- get(tbname)
    if(is.data.frame(item)){
      raw[[tbname]] <- item
    }
  } 
  
  rm(raw_table_names, tbname, item, sourmash_reports, mobtyper_reports, mlst_reports, mge_reports, metadata, contig_reports, amrfinder_reports, combined_amr_report)
}

# Set a proper "genome" key in each df

# Now this is handled in annotation_scan.R [Aurélie - 12/07/2023]
# For each annotation output data, an appropriate genome column is added
# See 'date'_annotation_report.xlsx 

# 2025-07-02
# Adding a merge especially for amr results and combined amr results
cat("Merge equivalent columns in AMR annotation results.\n")
{ 
  # Protein identifier == Protein id
  # Gene symbol == Element symbol
  # Sequence name == Element name
  # Element type == Type
  # Element subtype == Subtype
  
  raw$amrfinder_reports$`Protein identifier` <- ifelse(!is.na(raw$amrfinder_reports$`Protein identifier`), raw$amrfinder_reports$`Protein identifier`, raw$amrfinder_reports$`Protein id`)
  raw$amrfinder_reports$`Gene symbol` <- ifelse(!is.na(raw$amrfinder_reports$`Gene symbol`), raw$amrfinder_reports$`Gene symbol`, raw$amrfinder_reports$`Element symbol`)
  raw$amrfinder_reports$`Sequence name` <- ifelse(!is.na(raw$amrfinder_reports$`Sequence name`), raw$amrfinder_reports$`Sequence name`, raw$amrfinder_reports$`Element name`)
  raw$amrfinder_reports$`Element type` <- ifelse(!is.na(raw$amrfinder_reports$`Element type`), raw$amrfinder_reports$`Element type`, raw$amrfinder_reports$Type)
  raw$amrfinder_reports$`Element subtype` <- ifelse(!is.na(raw$amrfinder_reports$`Element subtype`), raw$amrfinder_reports$`Element subtype`, raw$amrfinder_reports$Subtype)
  
  raw$combined_amr_report$`Protein identifier` <- ifelse(!is.na(raw$combined_amr_report$`Protein identifier`), raw$combined_amr_report$`Protein identifier`, raw$combined_amr_report$`Protein id`)
  raw$combined_amr_report$`Gene symbol` <- ifelse(!is.na(raw$combined_amr_report$`Gene symbol`), raw$combined_amr_report$`Gene symbol`, raw$combined_amr_report$`Element symbol`)
  raw$combined_amr_report$`Sequence name` <- ifelse(!is.na(raw$combined_amr_report$`Sequence name`), raw$combined_amr_report$`Sequence name`, raw$combined_amr_report$`Element name`)
  raw$combined_amr_report$`Element type` <- ifelse(!is.na(raw$combined_amr_report$`Element type`), raw$combined_amr_report$`Element type`, raw$combined_amr_report$Type)
  raw$combined_amr_report$`Element subtype` <- ifelse(!is.na(raw$combined_amr_report$`Element subtype`), raw$combined_amr_report$`Element subtype`, raw$combined_amr_report$Subtype)
  
}

###############################
#* Gather data for each isolate as a JSON-like list
#* 
#* 
#* 

cat("Gathering data for each isolate as JSON_like list.\n")
{
  studied_genomes <- sort(intersect(raw$metadata$genome, raw$amrfinder_reports$genome))
  reports <- rbindlist(sapply(studied_genomes, function(i) {
    sub <- lapply(raw, function(tb) tb[genome == i])
    
    # How many resistance genes on chromosome ?
    
    amrf <- sub$amrfinder_reports
    ctg <- sub$contig_reports 
    amrf <- merge(amrf, ctg, by.x = "Contig id", by.y = "contig_id", all.x = TRUE, sort = FALSE)
    
    data.table(
      genome = i,
      cluster_id = sub$metadata$CLUSTER_ID,
      glims = sub$metadata$GLIMS,
      sampling_date = ymd(sub$metadata$DATE_PRELEVEMENT),
      ipp = sub$metadata$IPP,
      type_prlvt = sub$metadata$TYPE_PRELEVEMENT,
      soc_species = sub$metadata$SPECIES,
      sourmash_status = sub$sourmash_reports$status,
      sourmash_genus = sub$sourmash_reports$genus,
      sourmash_species = sub$sourmash_reports$species,
      mlst_species = sub$mlst_reports$mlst_species,
      sequence_type = sub$mlst_reports$sequence_type,
      args_chromosome = paste(sort(unique(amrf[`Element type` == "AMR" & molecule_type == "chromosome"]$`Gene symbol`)), collapse = ","),
      args_plasmid = paste(sort(unique(amrf[`Element type` == "AMR" & molecule_type == "plasmid"]$`Gene symbol`)), collapse = ","),
      plasmid_primary_clusters = paste(sort(unique(ctg[molecule_type == "plasmid"]$primary_cluster_id)), collapse = ","),
      plasmid_secondary_clusters = paste(sort(unique(ctg[molecule_type == "plasmid"]$secondary_cluster_id)), collapse = ","),
      plasmid_ID = paste(sort(unique(paste(ctg[molecule_type == "plasmid"]$primary_cluster_id, ctg[molecule_type == "plasmid"]$secondary_cluster_id, sep = "-"))), collapse = ","),
      args_primary_clusters = paste(amrf[`Element type` == "AMR" & molecule_type == "plasmid"]$primary_cluster_id, amrf[`Element type` == "AMR" & molecule_type == "plasmid"]$`Gene symbol`, sep = ":", collapse = ","),
      amr_classes = paste(sort(unique(str_to_title(amrf[`Element type` == "AMR"]$Class))), collapse = ","),
      betalactam_classes = paste(sort(unique(str_to_title(amrf[Class == "BETA-LACTAM"]$Subclass))), collapse = ",")
    )
  }, simplify = FALSE))
  
  rm(raw, studied_genomes)
}

cat("Saving to output file.\n")
{
  today_date <- Sys.Date()
  save(reports, file = paste0(today_date, "_", project, "_annotation_report.Rdata")) 
}

# Some cleaning
rm(reports, output_dir, project, today_date, oapply)

