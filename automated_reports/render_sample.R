####################
#### RENDERING #####
####################
## v.1.0
## 22.08.2023
##
rm(list = ls())
library(data.table)

## PARAMETERS ##
sample <- "023113373501-01"
illumina_assembly <- TRUE
outdir <- "~/Projets/4.Studies/2.Rapports/17.Epi-843"
path_scan <- "~/Projets/4.Studies/2.Rapports/17.Epi-843/2023-11-30_annotation_scan.Rdata"
path_reports <- "~/Projets/4.Studies/2.Rapports/17.Epi-843/2023-11-30_Epitrack_annotation_report.Rdata"
# outdir <- "../Epitrack_clusters/Samples"
# path_scan <- "../Script_Tests/2023-11-30_annotation_scan.Rdata"
# path_reports <- "../Script_Tests/2023-11-30_Epitrack_annotation_report.Rdata"
force <- TRUE

# Step1- Rendering sample report
if(illumina_assembly){
  out <- paste0(outdir, "/", sample, "_illumina_report.pdf")
} else {
  out <- paste0(outdir, "/", sample, "_nano_illumina_report.pdf")
}

if(!file.exists(out) | force == TRUE){
  rmarkdown::render(
    '~/Projets/2.Coding_projects/3.R/1.Automated_reports/Report/Non-interactive_sample_report.Rmd', 
    output_file = out,
    params = list(sample_glims = sample, 
                  illumina_only = illumina_assembly,
                  path_to_reports = path_reports,
                  path_to_scan = path_scan)
  )
}

