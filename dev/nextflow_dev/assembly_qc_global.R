########################
## ASSEMBLY QC GLOBAL ##
## 27-09-2023 ##########
########################

# Goal: Script appliable to Nano-Illumina assemblies and Illumina-only assemblies at the same time

rm(list = ls())
library(data.table)
library(stringr)
setwd("./genomes")

cat("Loading samples...\n")
# Nano-Illumina samples
nano_illumina_samples <- dir(".", "*polished.fasta$", recursive = TRUE)
nano_illumina_samples <- nano_illumina_samples[!grepl("/polish/", nano_illumina_samples)]
nano_illumina_samples <- str_extract(nano_illumina_samples, "^[^/]+")
cat(sprintf("There are %d Nano-Illumina samples. \n", length(nano_illumina_samples)))
# Illumina-only samples
illumina_samples <- list.dirs(".", recursive = FALSE, full.names = FALSE)
illumina_samples <- illumina_samples[!(illumina_samples %in% nano_illumina_samples)]
cat(sprintf("And there are %d Illumina-only samples.\n", length(illumina_samples)))

nano_illumina_samples
illumina_samples

cat("Assessing Illumina-only data.\n")
{
if(length(illumina_samples)!= 0){
  cat("    Scanning MLST data.\n")
  # Species info (from MLST)
  tmp <- list()
  for(s in illumina_samples) {
    fname <- sprintf("./%s/mlst/mlst.tsv", s)
    if(file.exists(fname))
    #x <- scan(fname, what="character", sep="\n")
    x <- read.table(fname, sep = "\t")
    species <- x[2]
    colnames(species) <- "species"
    tmp[[s]] <- data.table(id = s, species)
  }
  sample_species <- rbindlist(tmp)
  cat("    Scanning QUAST report for N50, number of contigs and length of largest contig.\n")
  # Assembly infos
  tmp <- list()
  for(s in illumina_samples) {
    fname <- sprintf("./%s/quast/transposed_report.tsv", s)
    tmp[[s]] <- cbind(data.table(id = s), fread(fname, sep = "\t"))
  }
  assemblies <- rbindlist(tmp)
  # We want no. of contigs, largest contig, coverage of largest contig
  assemblies_summary <- assemblies[ , .(
    N50 = N50,
    n_contig = `# contigs`,
    length_largest = `Largest contig`
  ) , by = id]
  cat("    Scanning Circlator report to get largest contig coverage.\n")
  # Assembly info step2
  # NB: Node_1 = largest contig
  tmp <- list()
  for(s in illumina_samples) {
    fname <- sprintf("./%s/circlator/%s_realigned.log", s, s)
    x <- read.table(fname, header = TRUE, sep = "\t")
    node_1 <- x$id[grepl("NODE_1_", x$id)]
    node_1_val <- str_split(node_1, "_")[[1]]
    cov_largest <- node_1_val[length(node_1_val)]
    tmp[[s]] <- data.table(id = s, cov_largest)
  }
  assemblies_2 <- rbindlist(tmp)
  cat("    Scanning Trimmomatic report to get number of reads.\n")
  # Reads infos (nb of reads)
  tmp <- list()
  for(s in illumina_samples) {
    fname <- sprintf("./%s/trimmomatic/trimmomatic.err", s)
    x <- scan(fname, what="character", sep="\n", quiet = TRUE)
    # print(str_extract(x[10], "[0-9]+"))
    n_reads <- as.double(str_extract(x[10], "[0-9]+"))
    tmp[[s]] <- data.table(id = s, n_reads)
    
  }
  reads <- rbindlist(tmp)
  cat("    Scanning SPAdes report to get average assembly coverage.\n")
  # Mean cov
  tmp <- list()
  for(s in illumina_samples) {
    fname <- sprintf("./%s/spades/spades_1.log", s)
    x <- scan(fname, what="character", sep="\n", quiet = TRUE)
    line_of_interest <- x[grepl("Estimated coverage", x) == TRUE][1]
    line_splitted <- unlist(str_split(line_of_interest, " "))
    mean_cov <- as.numeric(line_splitted[length(line_splitted)])
    tmp[[s]] <- data.table(id = s, mean_cov)
  }
  mean_coverage <- rbindlist(tmp)
  cat("    Merging results and creating output files.\n")
  output <- merge(sample_species, assemblies_summary)
  output <- merge(output, assemblies_2)
  output <- merge(output, reads)
  output <- merge(output, mean_coverage)
  # NB: n_contig = contig >500bp
  # Illumina: At least 500 000 reads, and maximum number of contigs = 150 and N50 >= 100 000
  output[ , QC_illumina_pass := FALSE ][ (species=="efaecium" & N50 >=30000 & n_contig <= 220 & n_reads >= 500000) | (N50 >= 100000 & n_contig <=150 & n_reads >= 500000), QC_illumina_pass := TRUE]
  
  # Writing to a file
  fwrite(output, file = "assembly_summary_illumina_only.tsv", sep = "\t")
  today_date <- Sys.Date()
  # Strains that pass QC
  cat(paste(output[QC_illumina_pass == TRUE]$id, collapse = "\n"), file = paste(today_date,"illumina_only_QC_passed.txt"))
  # Strains that fail Illumina QC
  cat(paste(output[QC_illumina_pass == FALSE]$id, collapse = "\n"), file = paste(today_date,"illumina_only_QC_illumina_failed.txt"))
  
  cat("    Some cleaning...\n")
  rm(assemblies, assemblies_2, assemblies_summary, output, reads, sample_species, species, tmp, cov_largest,
     fname, n_reads, node_1, node_1_val, s, today_date, x, mean_cov, line_splitted, line_of_interest, mean_coverage)
}
}

cat("Assessing Nano-Illumina data.\n")
{
if(length(nano_illumina_samples) != 0){
  cat("    Scanning MLST data.\n")
  # Species info (from MLST)
  tmp <- list()
  for(s in nano_illumina_samples) {
    fname <- sprintf("./%s/mlst/mlst.tsv", s)
    x <- read.table(fname, sep = "\t")
    species <- x[2]
    colnames(species) <- "species"
    tmp[[s]] <- data.table(id = s, species)
  }
  sample_species <- rbindlist(tmp)
  cat("    Scanning Flye report to get number of contigs, number of circular contigs, largest contig length and coverage.\n")
  # Assembly infos
  tmp <- list()
  for(s in nano_illumina_samples) {
    fname <- sprintf("./%s/flye/assembly_info.txt", s)
    tmp[[s]] <- cbind(data.table(id = s), fread(fname, sep = "\t"))
  }
  assemblies <- rbindlist(tmp)
  # We want no. of contigs, no. of circular contigs, largest contig, covergae of largest contig
  assemblies_summary <- assemblies[ , .(
    n_contig = .N,
    n_circular = sum(circ. == "Y"),
    length_largest = max(length),
    cov_largest = cov.[which.max(length)]
  ) , by = id]
  cat("    Scanning Bowtie2 mapping output to get overall and exact alignment rate (Illumina over Nanopore).\n")
  # Remapping info
  tmp <- list()
  for(s in nano_illumina_samples) {
    fname <- sprintf("./%s/polish/bowtie2.map.log", s)
    x <- scan(fname, what="character", sep="\n", quiet = TRUE)
    # print(x)
    illumina_reads <- as.double(str_extract(x[1], "^[0-9]+"))
    overall_alignment <- as.double(str_extract(x[15], "^[0-9.]+")) / 100
    # Idiom (?<=foo) is a non-matching group, not returned by extract. Similar to the (?:foo) non-capturing group
    exact_alignment <- as.double(str_extract(x[4], "(?<=\\()[0-9.]+")) /100
    tmp[[s]] <- data.table(id = s, illumina_reads, overall_alignment, exact_alignment)
  }
  alignments <- rbindlist(tmp)
  cat("    Scanning NanoPlot report for statistics on Nanopore sequencing data.\n")
  # Remapping info
  tmp <- list()
  for(s in nano_illumina_samples) {
    fname <- sprintf("./%s/nanoplot/NanoStats.txt", s)
    if(file.exists(fname)){
      x <- read.table(fname, sep = "\t", header = TRUE)
      nano_reads <- x[x$Metrics=="number_of_reads",2]
      nano_bases <- x[x$Metrics=="number_of_bases",2]
      nano_median_read_length <- x[x$Metrics=="median_read_length",2]
      nano_mean_read_length <- x[x$Metrics=="mean_read_length",2]
    } else {
      nano_reads <- NA
      nano_bases <- NA
      nano_median_read_length <- NA
      nano_mean_read_length <- NA
    }
    tmp[[s]] <- data.table(id = s, nano_reads, nano_bases, nano_median_read_length, nano_mean_read_length)
  }
  nano_stats <- rbindlist(tmp)
  cat("    Scanning Flye.log report for N50/N90.\n")
  # Remapping info
  tmp <- list()
  for(s in nano_illumina_samples) {
    fname <- sprintf("./%s/flye/flye.log", s)
    if(file.exists(fname)){
      x <- scan(fname, what="character", sep="\n", quiet = TRUE)
      line_of_interest <- x[grepl("Reads N50/N90", x) == TRUE][1]
      line_splitted <- unlist(str_split(line_of_interest, " "))
      N90 <- as.numeric(line_splitted[length(line_splitted)])
      N50 <- as.numeric(line_splitted[length(line_splitted)-2])
    } else {
      N50 <- NA
      N90 <- NA
    }
    tmp[[s]] <- data.table(id = s, N50, N90)
  }
  n50_stats <- rbindlist(tmp)
  cat("    Scanning Pilon report for polishing coverage and number of corrections.\n")  
  # Polishing info
  tmp <- list()
  for(s in nano_illumina_samples) {
    # Read log to extract coverage
    fname <- sprintf("./%s/polish/pilon.log", s)
    x <- scan(fname, what="character", sep="\n", quiet = TRUE)
    pilon_coverage <- as.double(str_extract(x[length(x)], "[0-9]+$"))
    # Read .changes file to extract polishing changes
    fname2 <- sprintf("./%s/polish/%s_polished.changes", s, s)
    x2 <- fread(fname2, header = FALSE)
    if(nrow(x2) > 0) {
      x2 <- x2[, .(V3, V4)]
      # Exclude long changes
      x2 <- x2[ str_length(V3) < 4 & str_length(V4) < 4]
      x2[ , .N, .(V3, V4)][order(-N)]
    } else print(s)
    pilon_corrections <- nrow(x2)
    tmp[[s]] <- data.table(id = s, pilon_coverage, pilon_corrections)
  }
  pilons <- rbindlist(tmp)
  cat("    Merging results and creating output files.\n")
  output <- merge(sample_species, assemblies_summary)
  output <- merge(output, nano_stats)
  output <- merge(output, n50_stats)
  output <- merge(output, alignments)
  output <- merge(output, pilons)

  # Apply QC checks then generate a list of ready nano_illumina_samples for the next step
    # Nanopore: At most 20 contigs, with coverage >= 20 on largest contig
  output[ , QC_nano_pass := FALSE ][ n_contig <= 20, QC_nano_pass := TRUE]
  
  # Illumina: At least 500 000 reads, 0.9 of which align on the assembly with a mean coverage of 30
  output[ , QC_illumina_pass := FALSE ][ pilon_coverage >= 30 & pilon_corrections >= 1e2 & pilon_corrections <= 1e4, QC_illumina_pass := TRUE]
  fwrite(output, file = "assembly_summary_nano_illumina.tsv", sep = "\t")
  # Excel-compatible output with forced text (avoid treating some fields as date, see https://stackoverflow.com/questions/165042/stop-excel-from-automatically-converting-certain-text-values-to-dates)
  output_xlsx <- data.table(output)
  output_xlsx[, id := sprintf('="%s"', id)]
  fwrite(output_xlsx, file = "assembly_summary_nano_illumina_excel.tsv", sep = "\t")
  currentDate <- Sys.Date()
  # Strains that pass both QC
  cat(paste(output[QC_nano_pass == TRUE & QC_illumina_pass == TRUE]$id, collapse = "\n"), file = paste(currentDate, "nano_illumina_QC_passed.txt"))
  # Strains that fail Nanopore QC
  cat(paste(output[QC_nano_pass == FALSE]$id, collapse = "\n"), file = paste(currentDate, "nano_illumina_QC_nanopore_failed.txt"))
  # Strains that fail Illumina QC
  cat(paste(output[QC_illumina_pass == FALSE]$id, collapse = "\n"), file = paste(currentDate, "nano_illumina_QC_illumina_failed.txt"))
  
  cat("    Some cleaning...\n")
  rm(alignments, assemblies, assemblies_summary, output, output_xlsx, pilons, sample_species, species, tmp, x2,
     currentDate, exact_alignment, fname, fname2, illumina_reads, nano_bases, nano_reads, overall_alignment, 
     pilon_corrections, pilon_coverage, s, x, nano_stats, n50_stats, line_splitted, line_of_interest)
}
}

rm(illumina_samples, nano_illumina_samples)
