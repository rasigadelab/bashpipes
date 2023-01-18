#######################################
# DNA READS PREPARATION ROUTINE

# Input: read files from possibly repeated Illumina and ONT runs, with a correspondence table linking file paths with sample

# Output: per-sample directory containing merged R1, R2 and ONT fastq.gz files ready for processing

# Clear R environment
rm(list = ls())
# Import R libraries
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
libraries <- c("data.table", "stringr", "argparse")
for (lib in libraries){usePackage(lib)}

# Add input parameters to script
parser <- ArgumentParser(description='Bactopia - Preparation of Data : merge of reads', add_help=TRUE)
parser$add_argument("-d", dest="path_to_data", required=TRUE, help="Path to the folder containing Illumina and Nanopore data", nargs='+')
args <- parser$parse_args()

# Set working directory
setwd(args$path_to_data)
# Create a directory for analysis
system("mkdir -p genomes")

# Remark: keep wsl shell open for faster execution (if wsl not open, each wsl command launches then shuts down the subsystem)
 
# Change Windows's CRLF format for Linux's LF format (mandatory for each generated file)
dos2unix <- function(fname) system(sprintf("dos2unix %s", fname))
# Run a bash command
bash <- function(cmd) system(cmd)
# Give exec permission for a file
perm <- function(sh) bash(sprintf("chmod u+x %s", sh))
# Shorthand: run a shell script without arguments
runsh <- function(sh) bash(sprintf("./%s", sh))

# List of read filenames and paths including read type R1, R2, ONT, and sample names
smplist <- data.table(read.table("Files_location.tsv", header = TRUE, sep = "\t"))

############ SANITY CHECKS

# Check file content against available files

x <- dir(".", pattern = "*.fastq.gz", recursive = TRUE)
sample_paths <- gsub("^([^/]*/[^/]*)/", "", smplist$full_path)
available_files <- "available_files.txt"
files <- c("Files that are in Files_location.tsv but not in Illumina or Nanopore folder.")
for(sample in setdiff(sample_paths, x)) files <- c(files, sample)
files <- c(files, "Files that are in Illumina or Nanopore folder, but not in Files_location.tsv.")
for(sample in setdiff(x, sample_paths)) files <- c(files, sample)
cat(paste(files, collapse = "\n"), file = "files_to_check.txt")

# ISSUE some files have invalid extension .fastq (2).gz
# EXAMPLE Illumina/ResisTrack_14/14-10_S4_R1.fastq (2).gz
# INFO later date and different file size compared to valid fastq.gz file in same directory
# CAUSE different runs pushed in same folder ?
# PREVENTION don't modify file extension; don't copy/paste without manual check; always append run ID with sample ID in filename
 
# FIX Ignore them for now

############ READ PREPARATION

# Generate merger commands

# Identify Illumina-only, Nanopore-only and hybrid samples

smplist[, has_illumina := any(type %in% c("R1", "R2")), by = Sample]
smplist[, has_nanopore := any(type %in% c("ONT")), by = Sample]

# Check filenames valid (no space)
x <- which(grepl(" ", smplist$full_path))
# smplist[x]

# Exclude them for now on (if length(x)==0, smplist[-x] empty the table)
files <- c(files,"Files that contain spaces in their name.")
if (length(x)!=0){
  for(sample in smplist[x]) files <- c(files, sample)
    smplist <- smplist[-x]
}
cat(paste(files, collapse = "\n"), file = "files_to_check.txt")
 
# Make directories
# TODO stop with warning if directories exist ?

workdir <- "genomes"

samples <- sort(unique(smplist$Sample))

mkdir.sh <- "makedirectories.sh"
mkdir.cmd <- "#!/bin/bash"
mkdir.cmd <- c(mkdir.cmd, "cd genomes")
for(sample in samples) mkdir.cmd <- c(mkdir.cmd, sprintf("mkdir -p %s", sample))
cat(paste(mkdir.cmd, collapse = "\n"), file = mkdir.sh)
dos2unix(mkdir.sh)
perm(mkdir.sh)
runsh(mkdir.sh)

# Prepare merger commands; no need to separate merger and copy,
# they take the same time for large files, see https://unix.stackexchange.com/questions/305594/copying-faster-than-cp

# List files for each sample and type
t <- smplist[, .(files = paste(gsub("^([^/]*/[^/]*)/", "", full_path), collapse = " ")), by = .(Sample, type)]

# t
# Prepare cat command
t[, cmd := sprintf("cat %s > %s/%s/%s_%s.fastq.gz", files, workdir, Sample, Sample, type)]

merge.sh <- "merge.sh"
merge.cmd <- c("#!/bin/bash", t$cmd)
cat(paste(merge.cmd, collapse = "\n"), file = merge.sh)
dos2unix(merge.sh)
perm(merge.sh)
runsh(merge.sh)

# Alternative, try massive parallelization ? No real benefit, looks like a disk bound process

# merge.sh <- "merge.sh"
# merge.cmd <- c("#!/bin/bash", sprintf("%s &", t$cmd))
# cat(paste(merge.cmd, collapse = "\n"), file = merge.sh)
# dos2unix(merge.sh)
# runsh(merge.sh)

####### TMP

# Prepare batches of 30 samples to be run sequentially, avoiding RAM overflow

batches <- split(samples, ceiling(seq_along(samples)/30))

for(batch in names(batches)) {
 fname <- sprintf("samples_%s.txt", batch)
 cat(paste(batches[[batch]], collapse = "\n"), file = fname)
 dos2unix(fname)
}

# Some info on how to chain commands with process IDs

# https://unix.stackexchange.com/questions/68693/how-to-know-if-a-background-job-is-finished

# some more info to give name to script processes,
# https://stackoverflow.com/questions/3075682/how-to-set-the-process-name-of-a-shell-script

# and kill them by PID using ps command (ps -U jpr)

# Monitor RAM usage, eg, free --giga --seconds 5 > ram.log