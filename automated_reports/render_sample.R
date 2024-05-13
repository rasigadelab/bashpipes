####################
#### RENDERING #####
####################
## v.1.0
## 22.08.2023
##
rm(list = ls())
library(data.table)

## PARAMETERS ##
sample <- "" # sample id
illumina_assembly <- TRUE # FALSE if hybrid assembly has been made on sample, else TRUE
outdir <- "out" 
path_scan <- "2023-11-30_annotation_scan.Rdata" # Path to Annotation_scan file
path_reports <- "2023-11-30_Epitrack_annotation_report.Rdata" # Path to Annotation_report file
force <- TRUE # Reprint the report even if an old one already exist

# Step1- Rendering sample report
if(illumina_assembly){
  out <- paste0(outdir, "/", sample, "_illumina_report.pdf")
} else {
  out <- paste0(outdir, "/", sample, "_nano_illumina_report.pdf")
}

if(!file.exists(out) | force == TRUE){
  rmarkdown::render(
    'Non-interactive_sample_report.Rmd', 
    output_file = out,
    params = list(sample_glims = sample, 
                  illumina_only = illumina_assembly,
                  path_to_reports = path_reports,
                  path_to_scan = path_scan)
  )
}

