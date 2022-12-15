#######################################
# DNA READS PREPARATION ROUTINE

# Input: read files from possibly repeated Illumina and ONT runs, with a correspondence table linking file paths with sample

# Output: per-sample directory containing merged R1, R2 and ONT fastq.gz files ready for processing


library(data.table)
library(readxl)
library(stringr)

# Remark: keep wsl shell open for faster execution (if wsl not open, each wsl command launches then shuts down the subsystem)

# Change Windows's CRLF format for Linux's LF format (mandatory for each generated file)
dos2unix <- function(fname) shell(sprintf("wsl.exe -d Ubuntu dos2unix %s", fname))
# Run a wsl command
wsl <- function(cmd) shell(sprintf("wsl.exe -d Ubuntu %s", cmd))
# Shorthand: run a shell script without arguments
runsh <- function(sh) wsl(sprintf("./%s", sh))

# List of read filenames and paths including read type R1, R2, ONT, and sample names
smplist <- data.table(read_excel("Files_location.xlsx"))

############ SANITY CHECKS

# Check file content against available files

x <- dir(".", pattern = "*.fastq.gz", recursive = TRUE)
setdiff(smplist$full_path, x)
setdiff(x, smplist$full_path)

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
smplist[x]

# Exclude them for now on
smplist <- smplist[-x]


# Make directories
# TODO stop with warning if directories exist ?

workdir <- "genomes"

samples <- sort(unique(smplist$Sample))

mkdir.sh <- "makedirectories.sh"
mkdir.cmd <- "#!/bin/bash"
mkdir.cmd <- c(mkdir.cmd, "cd genomes")
for(sample in samples) mkdir.cmd <- c(mkdir.cmd, sprintf("mkdir %s", sample))
cat(paste(mkdir.cmd, collapse = "\n"), file = mkdir.sh)
dos2unix(mkdir.sh)
runsh(mkdir.sh)

# Prepare merger commands; no need to separate merger and copy,
# they take the same time for large files, see https://unix.stackexchange.com/questions/305594/copying-faster-than-cp

# List files for each sample and type
t <- smplist[, .(files = paste(full_path, collapse = " ")), by = .(Sample, type)]

t
# Prepare cat command
t[, cmd := sprintf("cat %s > %s/%s/%s_%s.fastq.gz", files, workdir, Sample, Sample, type)]

merge.sh <- "merge.sh"
merge.cmd <- c("#!/bin/bash", t$cmd)
cat(paste(merge.cmd, collapse = "\n"), file = merge.sh)
dos2unix(merge.sh)
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




