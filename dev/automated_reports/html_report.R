#*******************************************************************
# Title: html_report.R
# Description: Renders all HTML reports for isolates listed in routine file.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Aurélie Fischer

rm(list = ls())
library(readxl)

# Load necessary data
routine_list <- read.table("./2026-03-23_all_epitrack_to_assemble.txt", sep = ";", header = FALSE)

# Load metadata
metadata <- read_excel("X:/ResisTrack/Suivi/Metadata_eq_rasigade.xlsx", sheet = 1)

# Obtain unique filtered IDs
list_glims <- unique(metadata[metadata$SAMPLE_ID %in% routine_list$V1,]$GLIMS)
#list_glims <- unique(reports$glims[grepl("Epi", reports$genome) | grepl("PA", reports$genome)])

# Change exceptions that have Nanopore sequencing but do have Illumina-only report
exceptions_nano <- c("023041259701-01", "022186972001-01", "022178399701-03",
                     "022153065501-02", "022143495501-01", "022171124807-01")
metadata[metadata$GLIMS %in% exceptions_nano,]$NANOPORE_STATUT <- "NR"

# Create log
log_file <- file("reports_creation.log", open = "at")
sink(log_file, type = "message")


# Loop for report generation
path_reports <- "data/2026-03-26_Epitrack_annotation_report.Rdata"
path_scan <- "data/2026-03-26_annotation_scan.Rdata"
path_qc <- "qc"

for (idglims in list_glims) {
  print(idglims)
  # Check Nanopore status
  nanopore <- metadata[metadata$GLIMS==idglims & !is.na(metadata$GLIMS),]$NANOPORE_STATUT
  if(TRUE){ #!file.exists(sprintf("reports/%s.html", idglims))){
    cat("Génération du rapport pour :", idglims, "\n", append = TRUE, file = log_file)
    if(!(is.na(nanopore)) && nanopore=="FAIT"){
      # Handle errors with tryCatch
      tryCatch({
        rmarkdown::render(
          input = "sample_report.Rmd",
          output_file = sprintf("reports/%s.html", idglims),
          params = list(
            sample_glims = idglims,
            illumina_only = FALSE,
            path_to_reports = path_reports,
            path_to_scan = path_scan,
            path_to_qc = path_qc
          )
        )
        cat("Rapport généré avec succès pour :", idglims, "\n", append = TRUE, file = log_file)
      }, error = function(e) {
        cat("Erreur lors de la génération du rapport pour :", idglims, "\n", append = TRUE, file = log_file)
        cat("Message d'erreur :", e$message, "\n", append = TRUE, file = log_file)
      })
    } else {
      # Handles errors with tryCatch
      tryCatch({
        rmarkdown::render(
          input = "sample_report.Rmd",
          output_file = sprintf("reports/%s.html", idglims),
          params = list(
            sample_glims = idglims,
            illumina_only = TRUE,
            path_to_reports = path_reports,
            path_to_scan = path_scan,
            path_to_qc = path_qc
          )
        )
        cat("Rapport généré avec succès pour :", idglims, "\n", append = TRUE, file = log_file)
      }, error = function(e) {
        cat("Erreur lors de la génération du rapport pour :", idglims, "\n", append = TRUE, file = log_file)
        cat("Message d'erreur :", e$message, "\n", append = TRUE, file = log_file)
      })
    }
  }
}