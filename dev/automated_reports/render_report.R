####################
#### RENDERING #####
####################
## v.1.0
## 22.08.2023
##
rm(list = ls())
library(data.table)

## WHAT PART SHOULD BE RUN ##
scanning <- FALSE
reporting <- FALSE
clustering <- FALSE

## PRELIMINARY STEPS ##
cat("Scanning annotation results + clustering samples. [OPTIONAL]\n")
{
  setwd("../Script_Tests")
  if(scanning==TRUE){
    #Step0- Loading annotations in .Rdata
    do.call(file.remove, list(grep("*.Rdata",list.files(".", full.names = TRUE), value = TRUE)))
    do.call(file.remove, list(grep("*_annotation_report.xlsx",list.files(".", full.names = TRUE), value = TRUE)))
    source("annotation_scan.R")
  }
  if(reporting==TRUE){
    source("annotation_report.R")
  }
  if(clustering==TRUE){
    source("clustering.R")
  }
  setwd("../Report")
}


## PARAMETERS ##
cluster <- ""
studied_species <- ""
env <- "" # could be either an IPP number or "" if there is no environmental sample
outdir <- paste0("../Epitrack_clusters/Epitrack_cluster_", cluster)
today_date <- Sys.Date()
path_scan <-  paste0(today_date, "_annotation_scan.Rdata")
path_reports <- paste0(today_date, "_Epitrack_annotation_report.Rdata")
path_clone <- paste0(today_date, "_clone_reports.Rdata")
path_qc <- ""
## WHAT PART SHOULD BE RUN ##
global_report <- TRUE
clonal_report <- TRUE
samples_reports <- FALSE

cat("Create output directory if it does not exist.\n")
{
  if (!dir.exists(outdir)) {dir.create(outdir)}
}

cat("Render global cluster report, including multi-species. [FOR INTERNAL USE]\n")
{
  if(global_report==TRUE){
    # Step1- Rendering global cluster report
    out <- paste0(outdir, "/cluster_",cluster, "_report.pdf")
    rmarkdown::render(
      'Non-interactive_detailed_cluster_report.Rmd',
      output_file = out,
      params = list(cluster_id = cluster,
                    path_to_reports = path_reports,
                    path_to_scan = path_scan)
    )
  }
}


cat("Rendering clustering report, per replicon and per species.\n")
{
  if(clonal_report==TRUE){
    for(sp in studied_species){
      out <- paste0(outdir, "/cluster_", cluster, "_", sp, "_chromosome.pdf")
      rmarkdown::render(
        'Non-interactive_quick_cluster_report.Rmd',
        output_file = out,
        params = list(cluster_id = cluster,
                      species = sp,
                      replicon = "Chromosome",
                      env_samples = env,
                      path_to_reports = path_clone,
                      path_to_scan = path_scan)
      )
    }
  }
}


cat("Rendering all samples report of cluster.\n")
{
  if(samples_reports==TRUE){
    # Step3- Rendering sample reports
    outdir <- "../Epitrack_clusters/Samples"
    # Get samples id
    load(path_reports)
    samples_cluster <- reports[cluster_id == cluster]$glims
    rm(reports)
    for(id in samples_cluster){
      #.rs.restartR()
      out <- paste0(outdir, "/", id, "_illumina_report.pdf")
      if(file.exists(out)){
        rmarkdown::render(
          'Non-interactive_sample_report.Rmd',
          output_file = out,
          params = list(sample_glims = id,
                        illumina_only = TRUE,
                        path_to_reports = path_reports,
                        path_to_scan = path_scan, 
                        path_to_qc = path_qc)
        )
      }
    }
  }
  
}


# Cleaning
rm(cluster, outdir, path_clone, 
   path_reports, path_scan, studied_species, today_date,
   clonal_report, clustering, global_report, reporting,
   samples_reports, scanning, env, id, samples_cluster)

