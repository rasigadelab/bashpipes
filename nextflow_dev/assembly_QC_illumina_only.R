# Load sample names
rm(list = ls())
library(data.table)
library(stringr)
setwd(".")

# Gather all samples with assembly, they have a *realigned.fasta file

samples <- dir(".", "*realigned.fasta$", recursive = TRUE)
samples <- samples[!grepl("/circlator/", samples)]
samples <- str_extract(samples, "^[^/]+")

# Species info (from MLST)
tmp <- list()
for(s in samples) {
  fname <- sprintf("./%s/mlst/mlst.tsv", s)
  #x <- scan(fname, what="character", sep="\n")
  x <- read.table(fname, sep = "\t")
  species <- x[2]
  colnames(species) <- "species"
  tmp[[s]] <- data.table(id = s, species)
}

tmp
sample_species <- rbindlist(tmp)
sample_species

# Assembly infos
tmp <- list()
for(s in samples) {
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

#Assembly info step2
#NB: Node_1 = largest contig
tmp <- list()
for(s in samples) {
  fname <- sprintf("./%s/circlator/%s_realigned.log", s, s)
  #x <- scan(fname, what="character", sep="\n")
  x <- read.table(fname, header = TRUE, sep = "\t")
  node_1 <- x$id[grepl("NODE_1_", x$id)]
  node_1_val <- str_split(node_1, "_")[[1]]
  cov_largest <- node_1_val[length(node_1_val)]
  tmp[[s]] <- data.table(id = s, cov_largest)
}

tmp
assemblies_2 <- rbindlist(tmp)
assemblies_2

# Reads infos (nb of reads)
tmp <- list()
for(s in samples) {
  fname <- sprintf("./%s/trimmomatic/trimmomatic.err", s)
  x <- scan(fname, what="character", sep="\n")
  # print(str_extract(x[10], "[0-9]+"))
  n_reads <- as.double(str_extract(x[10], "[0-9]+"))
  tmp[[s]] <- data.table(id = s, n_reads)
  
}

tmp
reads <- rbindlist(tmp)
reads

output <- merge(sample_species, assemblies_summary)
output <- merge(output, assemblies_2)
output <- merge(output, reads)


#NB: n_contig = contig >500bp
# Illumina: At least 500 000 reads, and maximum number of contigs = 100
output[ , QC_illumina_pass := FALSE ][ N50 >= 100000 & n_contig <=150 & n_reads >= 500000, QC_illumina_pass := TRUE]

# Writing to a file
fwrite(output, file = "assembly_summary.tsv", sep = "\t")
today_date <- Sys.Date()
# Strains that pass QC
cat(paste(output[QC_illumina_pass == TRUE]$id, collapse = "\n"), file = paste(today_date,"_QC_passed.txt"))
# Strains that fail Illumina QC
cat(paste(output[QC_illumina_pass == FALSE]$id, collapse = "\n"), file = paste(today_date,"_QC_illumina_failed.txt"))
