# Load sample names
library(data.table)
library(stringr)

# Gather all samples with assembly, they have a *polished.fasta file

samples <- dir(".", "*polished.fasta", recursive = TRUE)
samples <- samples[!grepl("/polish/", samples)]
samples <- str_extract(samples, "^[^/]+")

x <- character(0)
for(i in 1:10) x <- c(x, scan(sprintf("samples_%i.txt", i), what="character", sep="\n"))

setdiff(x, samples)

# samples <- read.csv("samples.txt", header = F)[[1]]

# Assembly infos
tmp <- list()
for(s in samples) {
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

# Remapping info


tmp <- list()
for(s in samples) {
  fname <- sprintf("./%s/polish/bowtie2.map.log", s)
  x <- scan(fname, what="character", sep="\n")
  # print(x)
  n_reads <- as.double(str_extract(x[1], "^[0-9]+"))
  overall_alignment <- as.double(str_extract(x[15], "^[0-9.]+")) / 100
  # Idiom (?<=foo) is a non-matching group, not returned by extract. Similar to the (?:foo) non-capturing group
  exact_alignment <- as.double(str_extract(x[4], "(?<=\\()[0-9.]+")) /100
  tmp[[s]] <- data.table(id = s, n_reads, overall_alignment, exact_alignment)
  
}

tmp
alignments <- rbindlist(tmp)
alignments

# Polishing info
tmp <- list()
for(s in samples) {
  # Read log to extract coverage
  fname <- sprintf("./%s/polish/pilon.log", s)
  x <- scan(fname, what="character", sep="\n")
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

tmp
pilons <- rbindlist(tmp)

output <- merge(assemblies_summary, alignments)
output <- merge(output, pilons)
output

boxplot(cov_largest ~ n_contig, output)
cor.test(output$cov_largest, output$overall_alignment, method = "s")
cor.test(output$cov_largest, output$exact_alignment, method = "s")

# Apply QC checks then generate a list of ready samples for the next step

# Nanopore: At most 20 contigs, with coverage >= 20 on largest contig
output[ , QC_nano_pass := FALSE ][ n_contig <= 20 & cov_largest >= 20, QC_nano_pass := TRUE]

# Illumina: At least 500 000 reads, 0.9 of which align on the assembly with a mean coverage of 30
# output[ , QC_illumina_pass := FALSE ][ n_reads  >= 5e5 & overall_alignment >= 0.9 & pilon_coverage >= 30, QC_illumina_pass := TRUE]
output[ , QC_illumina_pass := FALSE ][ pilon_coverage >= 30 & pilon_corrections >= 1e2 & pilon_corrections <= 1e4, QC_illumina_pass := TRUE]

# View(output)

mean(output$QC_nano_pass)
mean(output$QC_illumina_pass)

sum(output$QC_illumina_pass & output$QC_nano_pass)

hist(output[, .(ratio = exact_alignment / overall_alignment)]$ratio)
hist(output$exact_alignment)
hist(output$n_reads)

fwrite(output, file = "assembly_summary.tsv", sep = "\t")

# Excel-compatible output with forced text (avoid treating some fields as date, see https://stackoverflow.com/questions/165042/stop-excel-from-automatically-converting-certain-text-values-to-dates)
output_xlsx <- data.table(output)
output_xlsx[, id := sprintf('="%s"', id)]

fwrite(output_xlsx, file = "assembly_summary_excel.tsv", sep = "\t")


# Strains that pass both QC
cat(paste(output[QC_nano_pass == TRUE & QC_illumina_pass == TRUE]$id, collapse = "\n"), file = "221026_QC_passed.txt")
# Strains that fail Nanopore QC
cat(paste(output[QC_nano_pass == FALSE]$id, collapse = "\n"), file = "221026_QC_nanopore_failed.txt")
# Strains that fail Illumina QC
cat(paste(output[QC_illumina_pass == FALSE]$id, collapse = "\n"), file = "221026_QC_illumina_failed.txt")

