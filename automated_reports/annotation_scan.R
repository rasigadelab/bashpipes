#*******************************************************************
#*GENOME ANNOTATION SCAN
#*V0.2 2023-08-02
#*
#* Scans annotation files to feed the annotation database
#* 

rm(list = objects())

library(data.table)
library(stringr)
library(readr)
options(java.parameters = "-Xmx5000m") # Avoid heap overflow in Java, see https://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
library(openxlsx)

# Change Windows's CRLF format for Linux's LF format (mandatory for each generated file)
dos2unix <- function(fname) shell(sprintf("wsl.exe -d Ubuntu dos2unix %s", fname))
# Run a wsl command
wsl <- function(cmd) shell(sprintf("wsl.exe -d Ubuntu %s", cmd))
# Run a wsl command with embedded sprintf
wslf <- function(...) shell(sprintf("wsl.exe -d Ubuntu %s", sprintf(...)))
# Shorthand: run a shell script without arguments
runsh <- function(sh) wsl(sprintf("./%s", sh))


# PARAMETERS

# Where to write the results
output_dir <- getwd()

# Where to find the genomes
genome_dir <- "."

# Include an excel output ? (may be lengthy)
excel_output <- TRUE

# Do environment cleaning at the end ? 
clean_env <- TRUE


cat(sprintf("Entering genome directory: %s\n", genome_dir))

setwd(genome_dir)

# Analyse file names
fnames <- dir(recursive = TRUE, full.names = TRUE)
inames <- unique(str_match(fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"])

cat(sprintf("  Genome directory contains %i files in %i isolates.\n", length(fnames), length(inames)))

# TODO explicit key = genome

# Gather amrfinder and mob_recon informations then combine in data tables

# MOB-SUITE Contig_reports
cat("Scanning MOB_SUITE reports.\n")
{
  contig_report_fnames <- fnames[grepl("contig_report.txt", fnames)]
  contig_report_inames <- str_match(contig_report_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(contig_report_fnames) == length(contig_report_inames))
  names(contig_report_fnames) <- contig_report_inames
  cat(sprintf("  Found %i isolates with a MOB_SUITE CONTIG report.\n", length(contig_report_inames)))
   
  contig_reports <- list()
  for(iname in contig_report_inames) {
    contig_reports[[ iname ]] <- fread( contig_report_fnames[iname] )[ , sample_id := iname]
  }
  contig_reports <- rbindlist(contig_reports)[, key := paste(sample_id, contig_id, sep = "_")]
  setnames(contig_reports, "sample_id", "genome")
  
  rm(contig_report_fnames, contig_report_inames, iname)
}
cat("  Finished scanning MOB_SUITE CONTIG reports.\n")

# MOB-SUITE mge_reports
{
  mge_fnames <- fnames[grepl("mge.report.txt", fnames)]
  mge_inames <- str_match(mge_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(mge_fnames) == length(mge_inames))
  names(mge_fnames) <- mge_inames
  cat(sprintf("  Found %i isolates with a MOB_SUITE MGE report.\n", length(mge_inames)))
  
  mge_reports <- list()
  for(iname in mge_inames) {
    mge_reports[[ iname ]] <- fread( mge_fnames[iname] )[ , sample_id := iname]
  }
  mge_reports <- rbindlist(mge_reports)[, key := paste(sample_id, contig_id, sep = "_")]
  setnames(mge_reports, "sample_id", "genome")
  
  rm(mge_fnames, mge_inames, iname)
}
cat("  Finished scanning MOB_SUITE MGE reports.\n")

# MLST reports
cat("Scanning MLST reports.\n")
{
  mlst_fnames <- fnames[grepl("mlst.tsv", fnames)]
  mlst_inames <- str_match(mlst_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  names(mlst_fnames) <- mlst_inames
  cat(sprintf("  Found %i isolates with a MLST report.\n", length(mlst_inames)))
  
  mlst_reports <- list()
  for(iname in mlst_inames) {
    mlst_reports[[ iname ]] <- fread( mlst_fnames[iname], header = FALSE )[ , V1 := iname]
  }
  mlst_reports <- rbindlist(mlst_reports, fill = TRUE)#[, key := paste(sample_id, contig_id, sep = "_")]
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
  amrfinder_report_fnames <- fnames[grepl("amrfinder.tsv", fnames)]
  amrfinder_report_inames <- str_match(amrfinder_report_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(amrfinder_report_fnames) == length(amrfinder_report_inames))
  names(amrfinder_report_fnames) <- amrfinder_report_inames
  cat(sprintf("  Found %i isolates with an AMRFINDER+ report.\n", length(amrfinder_report_inames)))
  
  amrfinder_reports <- list()
  for(iname in amrfinder_report_inames) {
    # amrfinder_reports[[ iname ]] <- fread( amrfinder_report_fnames[iname] )[, genome := iname]
    amrfinder_reports[[ iname ]] <- fread( amrfinder_report_fnames[iname] )
    amrfinder_reports[[ iname ]] <- cbind(genome = iname, amrfinder_reports[[ iname ]])
  }
  amrfinder_reports <- rbindlist(amrfinder_reports)[, key := paste(genome, `Contig id`, sep = "_")]

  rm(amrfinder_report_fnames, amrfinder_report_inames, iname)
}
cat("  Finished scanning AMRFinder+ reports.\n")

# Sourmash reports
cat("Scanning SOURMASH reports.\n")
{
  sourmash_report_fnames <- fnames[grepl("sourmash.csv", fnames)]
  sourmash_report_inames <- str_match(sourmash_report_fnames, "^\\./(?<isolate>[^/]+)")[, "isolate"]
  stopifnot(length(sourmash_report_fnames) == length(sourmash_report_inames))
  names(sourmash_report_fnames) <- sourmash_report_inames
  cat(sprintf("  Found %i isolates with a SOURMASH report.\n", length(sourmash_report_inames)))
  
  sourmash_reports <- list()
  for(iname in sourmash_report_inames) {
    sourmash_reports[[ iname ]] <- fread( sourmash_report_fnames[iname] )[, ID := iname]
  }
  sourmash_reports <- rbindlist(sourmash_reports)
  setnames(sourmash_reports, "ID", "genome")
  
  rm(sourmash_report_fnames, sourmash_report_inames, iname)
}
cat("  Finished scanning SOURMASH reports.\n")

cat("Leaving genome directory.\n")
###END SCAN---
rm(fnames, genome_dir, inames)

cat(sprintf("Entering output directory: %s\n", output_dir))
setwd(output_dir)

# Merge reports
combined_amr_report <- merge(amrfinder_reports, contig_reports, by = c("key", "genome"), all.x = TRUE)
combined_amr_report <- merge(combined_amr_report, sourmash_reports, by = "genome", all.x = TRUE)
combined_amr_report <- merge(combined_amr_report, mlst_reports, by = "genome", all.x = TRUE)

# Just watching some results
# Sourmash detected genus
combined_amr_report[ !duplicated(genome) , .N, genus][ order(-N)]
# SOC species
combined_amr_report[ !duplicated(genome) , .N, species][ order(-N)]
# ARGs
combined_amr_report[ , .N, `Gene symbol`][ order(-N) ]
# Plasmid primary cluster ID
combined_amr_report[ !duplicated(paste(genome, primary_cluster_id)), .N, primary_cluster_id][ order(-N) ]

############################
# SAVE

cat("Saving R data format.\n")
today_date <- Sys.Date()
save(mge_reports, amrfinder_reports, sourmash_reports, contig_reports, mlst_reports, combined_amr_report, file = paste0(today_date, "_annotation_scan.Rdata"))

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
  rm(amrfinder_reports, combined_amr_report, contig_reports, data_to_save, mge_reports, mlst_reports, sourmash_reports, today_date, sheet_names)
}

rm(clean_env, excel_output, output_dir, dos2unix, runsh, wsl, wslf)

