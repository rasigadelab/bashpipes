########################
## ASSEMBLY QC GLOBAL ##
##     12-01-2025     ##
########################

# Goal: Script appliable to Nano-Illumina assemblies and Illumina-only assemblies at the same time

rm(list = ls())
library(data.table)
library(stringr)
library(readxl)
library(dplyr)
setwd("./genomes")
metadata_file <- ""

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
  cat("    Scanning QUAST report for N50, number of contigs, length of largest contig and assembly length.\n")
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
    length_largest = `Largest contig`,
    length_assembly = `Total length (>= 0 bp)`
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
    line_with_n_reads <- grep("Input Read", x)
    n_reads <- as.double(str_extract(x[line_with_n_reads], "[0-9]+"))
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
  cat("    Scanning Sourmash report to get species.\n")
  # Sourmash species
  tmp <- list()
  for(s in illumina_samples) {
    fname <- sprintf("./%s/sourmash/sourmash.csv", s)
    x <- read.table(fname, header = TRUE, sep=',')
    tmp[[s]] <- data.table(id = s, sourmash_species = x$species)
  }
  sourmash_species <- rbindlist(tmp)
  cat("    Merging results and creating output files.\n")
  output <- merge(sample_species, sourmash_species)
  output <- merge(output, assemblies_summary)
  output <- merge(output, assemblies_2)
  output <- merge(output, reads)
  output <- merge(output, mean_coverage)
  # NB: n_contig = contig >500bp
  # NEED TO WORK 
  # HERE 
  metadata <- data.table(read_excel(metadata_file, sheet=1))
  rm(metadata_file)
  
  # Merge both datatables to include MALDI species
  metadata <- metadata[,c("SAMPLE_ID", "SPECIES")]
  study_df <- merge(metadata, output, by.x = "SAMPLE_ID", by.y = "id", all.y = TRUE)
  colnames(study_df) <- c("genome", "MALDI_species", "MLST_species", "sourmash_species", 
                          "N50", "n_contig", "length_largest", "length_assembly",
                          "cov_largest", "n_reads", "mean_cov")
  
  ################################################################################################
  ## Two things to add : New QC pass and Contamination indication
  ################################################################################################
  ## This is based on QC criteria work end of November 2025
  ## QC criteria per species.
  ## Some genus with low number of samples were not included inside analysis. Old criteria will be kept for them.
  
  # Let's apply new QC thresholds
  # For all, n_reads>500 000
  # Rare genus
  rare_genus <- c("Achromobacter", "Alcaligenes", "Burkholderia", "Delftia", "Bacillus", "Bacteroides",
                  "Campylobacter", "Corynebacterium", "Turicella", "Hafnia", "Kluyvera", "Proteus",
                  "Raoultella", "Shigella", "Lactobacillus", "Moraxella", "Pasteurella", "Propionibacterium",
                  "Rhizobium", "Morganella", "Pseudocitrobacter")
  
  study_df <- study_df %>%
    # KLEBSIELLA
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Klebsiella oxytoca" & n_reads>500000 & n_contig<=147 & N50>=79000 & length_assembly>=5500000 & length_assembly<=6900000, TRUE, FALSE)) %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Klebsiella pneumoniae" & n_reads>500000 & n_contig<=162 & N50>=71000 & length_assembly>=5000000 & length_assembly<=6300000, TRUE, NEW_QC_illumina_pass)) %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Klebsiella variicola" & n_reads>500000 & n_contig<=109 & N50>=165000 & length_assembly>=6000000 & length_assembly<=6100000, TRUE, NEW_QC_illumina_pass)) %>%
    # RARE GENUS
    mutate(NEW_QC_illumina_pass=ifelse(str_extract(MALDI_species, "[A-Z][a-z]+") %in% rare_genus & n_reads>500000 & n_contig<=150 & N50>=100000, TRUE, NEW_QC_illumina_pass))  %>%
    # ACINETOBACTER
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Acinetobacter baumannii" & n_reads>500000 & n_contig<=122 & N50>=60000 & length_assembly>=3600000 & length_assembly<=4300000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Acinetobacter johnsonii" & n_reads>500000 & n_contig<=211 & N50>=34000 & length_assembly>=3400000 & length_assembly<=3800000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Acinetobacter nosocomialis" & n_reads>500000 & n_contig<=40 & N50>=197000 & length_assembly>=4000000 & length_assembly<=4100000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Acinetobacter ursingii" & n_reads>500000 & n_contig<=146 & N50>=122000 & length_assembly>=3300000 & length_assembly<=4300000, TRUE, NEW_QC_illumina_pass))  %>%
    # ESCHERICHIA
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Escherichia fergusonii" & n_reads>500000 , NA, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Escherichia coli" & n_reads>500000 & N50<50000, FALSE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Escherichia coli" & n_reads>500000 & n_contig<=245 & N50>=50000 & N50<=120000 & length_assembly>=4300000 & length_assembly<=5700000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Escherichia coli" & n_reads>500000 & n_contig<=151 & N50>120000 & length_assembly>=4500000 & length_assembly<=5600000, TRUE, NEW_QC_illumina_pass))  %>%
    # CITROBACTER
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter amalonaticus" & n_reads>500000 & n_contig<=210 & N50>=15000 & length_assembly>=4900000 & length_assembly<=5500000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter braakii" & n_reads>500000 & n_contig<=188 & N50>=53000 & length_assembly>=4800000 & length_assembly<=6200000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter farmeri" & n_reads>500000 & n_contig<=85 & N50>=18000 & length_assembly>=4800000 & length_assembly<=5800000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter freundii" & n_reads>500000 & n_contig<=183 & N50>=63000 & length_assembly>=5000000 & length_assembly<=6100000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter freundii / braakii" & n_reads>500000 & n_contig<=176 & N50>=100000 & length_assembly>=4700000 & length_assembly<=6200000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter freundii complex" & n_reads>500000 & n_contig<=217 & N50>=33000 & length_assembly>=5000000 & length_assembly<=6100000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter koseri" & n_reads>500000 & n_contig<=58 & N50>=200000 & length_assembly>=4100000 & length_assembly<=5500000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter murliniae" & n_reads>500000, NA, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter sedlakii" & n_reads>500000, NA, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter spp." & n_reads>500000 & n_contig<=138 & N50>=90000 & length_assembly>=5400000 & length_assembly<=6000000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter youngae" & n_reads>500000 & n_contig<=214 & N50>=100000 & length_assembly>=4600000 & length_assembly<=6300000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Citrobacter koseri-amalonaticus/farmeri" & n_reads>500000 , NA, NEW_QC_illumina_pass))  %>%
    # PSEUDOMONAS
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Pseudomonas aeruginosa" & n_reads>500000 & n_contig<=193 & N50>=50000 & length_assembly>=5600000 & length_assembly<=7600000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Pseudomonas citronellolis" & n_reads>500000 & n_contig<=390 & N50>=31000 & length_assembly>=6600000 & length_assembly<=7500000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Pseudomonas fluorescens" & n_reads>500000 & n_contig<=217 & N50>=29000 & length_assembly>=5500000 & length_assembly<=7100000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Pseudomonas putida" & n_reads>500000 & n_contig<=350 & N50>=30000 & length_assembly>=5000000 & length_assembly<=6900000, TRUE, NEW_QC_illumina_pass))  %>%
    # SERRATIA
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Serratia marcescens" & n_reads>500000 & n_contig<=250 & N50>=50000 & length_assembly>=4600000 & length_assembly<=6000000, TRUE, NEW_QC_illumina_pass))  %>%
    # ENTEROBACTER
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter aerogenes" & n_reads>500000 & n_contig<=92 & N50>=67000 & length_assembly>=4900000 & length_assembly<=5600000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter asburiae" & n_reads>500000 & n_contig<=131 & N50>=106000 & length_assembly>=4700000 & length_assembly<=5400000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter cloacae" & n_reads>500000 & n_contig<=137 & N50>=30000 & length_assembly>=4500000 & length_assembly<=5500000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter cloacae complex" & n_reads>500000 & n_contig<=140 & N50>=30000 & length_assembly>=4600000 & length_assembly<=5600000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter cloacae ssp. cloacae" & n_reads>500000 & n_contig<=104 & N50>=116000 & length_assembly>=5300000 & length_assembly<=5400000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter cloacae ssp. dissolvens" & n_reads>500000 & n_contig<=70 & N50>=188000 & length_assembly>=5000000 & length_assembly<=5100000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter hormaechei" & n_reads>500000 & n_contig<=67 & N50>=182000 & length_assembly>=4800000 & length_assembly<=5400000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter hormaechei ssp. hoffmannii" & n_reads>500000 & n_contig<=87 & N50>=116000 & length_assembly>=5100000 & length_assembly<=5200000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter kobei" & n_reads>500000 & n_contig<=119 & N50>=101000 & length_assembly>=5100000 & length_assembly<=5300000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterobacter roggenkampii" & n_reads>500000 & n_contig<=94 & N50>=136000 & length_assembly>=5300000 & length_assembly<=5400000, TRUE, NEW_QC_illumina_pass))  %>%
    # ENTEROCOCCUS
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterococcus faecalis" & n_reads>500000 & n_contig<=93 & N50>=100000 & length_assembly>=2600000 & length_assembly<=3500000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Enterococcus faecium" & n_reads>500000 & n_contig<=246 & N50>=25000 & length_assembly>=2800000 & length_assembly<=3100000, TRUE, NEW_QC_illumina_pass))  %>%
    # HAEMOPHILUS
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Haemophilus influenzae" & n_reads>500000 & n_contig<=53 & N50>=75000 & length_assembly>=1600000 & length_assembly<=2000000, TRUE, NEW_QC_illumina_pass))  %>%
    # STAPHYLOCOCCUS
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus aureus" & n_reads>500000 & n_contig<=76 & N50>=55000 & length_assembly>=2500000 & length_assembly<=3000000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus capitis" & n_reads>500000 & n_contig<=53 & N50>=71000 & length_assembly>=2500000 & length_assembly<=2700000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus caprae" & n_reads>500000 & n_contig<=44 & N50>=100000 & length_assembly>=2500000 & length_assembly<=2900000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus epidermidis" & n_reads>500000 & n_contig<=103 & N50>=51000 & length_assembly>=2200000 & length_assembly<=2800000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus haemolyticus" & n_reads>500000 & n_contig<=111 & N50>=51000 & length_assembly>=2400000 & length_assembly<=2600000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus hominis" & n_reads>500000 & n_contig<=81 & N50>=50000 & length_assembly>=2100000 & length_assembly<=2400000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus lugdunensis" & n_reads>500000 & n_contig<=62 & N50>=100000 & length_assembly>=1500000 & length_assembly<=2500000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus pettenkoferi" & n_reads>500000 & n_contig<=77 & N50>=102000 & length_assembly>=2400000 & length_assembly<=2600000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus saprophyticus" & n_reads>500000 & n_contig<=32 & N50>=151000 & length_assembly>=2600000 & length_assembly<=2700000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Staphylococcus simulans" & n_reads>500000 & n_contig<=69 & N50>=231000 & length_assembly>=2600000 & length_assembly<=2800000, TRUE, NEW_QC_illumina_pass))  %>%
    # STENOTROPHOMONAS
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Stenotrophomonas maltophilia" & n_reads>500000 & n_contig<=219 & N50>=25000 & length_assembly>=4300000 & length_assembly<=5200000, TRUE, NEW_QC_illumina_pass))  %>%
    # STREPTOCOCCUS
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Streptococcus anginosus" & n_reads>500000 & n_contig<=57 & N50>=100000 & length_assembly>=1600000 & length_assembly<=2400000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Streptococcus constellatus" & n_reads>500000 & n_contig<=102 & N50>=30000 & length_assembly>=1600000 & length_assembly<=2400000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Streptococcus dysgalactiae" & n_reads>500000 & n_contig<=64 & N50>=69000 & length_assembly>=2000000 & length_assembly<=2200000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Streptococcus pneumoniae" & n_reads>500000 & n_contig<=102 & N50>=36000 & length_assembly>=1800000 & length_assembly<=2200000, TRUE, NEW_QC_illumina_pass))  %>%
    mutate(NEW_QC_illumina_pass=ifelse(MALDI_species=="Streptococcus pyogenes" & n_reads>500000 & n_contig<=28 & N50>=98000 & length_assembly>=1600000 & length_assembly<=2000000, TRUE, NEW_QC_illumina_pass))
  
  ###########################################################################################################
  # Add a variable that indicates contamination
  # On checke d'abord si concordance entre Sourmash et MALDI
  # Si pas l'info de sourmash, alors on checke la concordance entre MLST et MALDI
  list_species <- sort(unique(study_df$MALDI_species))
  study_df <- study_df %>%
    # No MLST No Sourmash
    mutate(contaminated=ifelse(is.na(sourmash_species) & MLST_species=="-", "NO_PROOF", TRUE)) %>%
    # ACHROMOBACTER
    mutate(contaminated=ifelse(MALDI_species=="Achromobacter xylosidans" & sourmash_species %in% c("s__Achromobacter xylosoxidans", "s__Achromobacter sp002902905"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Achromobacter xylosidans" & is.na(sourmash_species) & MLST_species %in% c("achromobacter"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Achromobacter xylosidans/denitrificans" & sourmash_species %in% c("s__Achromobacter xylosoxidans", "s__Achromobacter insuavis_A"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Achromobacter xylosidans/denitrificans" & is.na(sourmash_species) & MLST_species %in% c("achromobacter"), FALSE, contaminated)) %>%
    # ACINETOBACTER
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter baumannii" & sourmash_species %in% c("s__Acinetobacter baumannii", "Acinetobacter baumannii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter baumannii" & is.na(sourmash_species) & MLST_species %in% c("abaumannii_2"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter johnsonii" & sourmash_species %in% c("s__Acinetobacter johnsonii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter johnsonii" & is.na(sourmash_species) & MLST_species %in% c("abaumannii_2"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter nosocomialis" & sourmash_species %in% c("s__Acinetobacter nosocomialis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter nosocomialis" & is.na(sourmash_species) & MLST_species %in% c("abaumannii_2"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter ursingii" & sourmash_species %in% c("s__Acinetobacter ursingii", "Acinetobacter ursingii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Acinetobacter ursingii" & is.na(sourmash_species) & MLST_species %in% c("abaumannii_2"), FALSE, contaminated)) %>%
    # ALCALIGENES
    mutate(contaminated=ifelse(MALDI_species=="Alcaligenes faecalis" & sourmash_species %in% c("Alcaligenes faecalis"), FALSE, contaminated)) %>%
    # BACILLUS
    mutate(contaminated=ifelse(MALDI_species=="Bacillus cereus" & sourmash_species %in% c("Bacillus cereus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Bacillus cereus" & is.na(sourmash_species) & MLST_species %in% c("bcereus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Bacillus thuringiensis" & sourmash_species %in% c("s__Bacillus_A thuringiensis_S", "s__Bacillus_A thuringiensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Bacillus thuringiensis" & is.na(sourmash_species) & MLST_species %in% c("bcereus"), FALSE, contaminated)) %>%
    # BACTEROIDES
    mutate(contaminated=ifelse(MALDI_species=="Bacteroides fragilis" & sourmash_species %in% c("Bacteroides fragilis"), FALSE, contaminated)) %>%
    # BURKHOLDERIA
    mutate(contaminated=ifelse(MALDI_species=="Burkholderia multivorans" & is.na(sourmash_species) & MLST_species %in% c("bcc"), FALSE, contaminated)) %>%
    # CAMPYLOBACTER
    mutate(contaminated=ifelse(MALDI_species=="Campylobacter jejunii" & sourmash_species %in% c("Campylobacter jejuni"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Campylobacter jejunii" & is.na(sourmash_species) & MLST_species %in% c("campylobacter"), FALSE, contaminated)) %>%
    # CITROBACTER 
    # normalement devrait se différencier : AMALONATICUS - FARMERI - KOSERI (chacun différent les uns des autres) de la clade FREUNDII COMPLEX [freundii-braakii-youngae-et-cie]
    # Mais MALDI discrimine mal les espèces des Citrobacter, apparently amalonaticus isok, mais farmeri pas très bien distingué
    # Après sur les outils génomiques : MLST ne discrimine pas du tout non plus, tout le monde est dans la case cfreundii
    # Sourmash se débrouille déjà mieux, il distingue pas très bien les C. youngae
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter amalonaticus" & sourmash_species %in% c("s__Citrobacter_A amalonaticus", "s__Citrobacter_A telavivensis", "s__Citrobacter freundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter amalonaticus" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter braakii" & sourmash_species %in% c("Citrobacter braakii", "Citrobacter freundii", "s__Citrobacter portucalensis", "s__Citrobacter braakii", "s__Citrobacter europaeus", "s__Citrobacter_A farmeri", "s__Citrobacter freundii", "s__Citrobacter murliniae", "Citrobacter sp. MGH103", "Citrobacter sp. A316"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter braakii" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter farmeri" & sourmash_species %in% c("s__Citrobacter_A farmeri"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter farmeri" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter freundii" & sourmash_species %in% c("Citrobacter freundii", "Citrobacter sp. MGH103", "s__Citrobacter braakii", "s__Citrobacter freundii", "s__Citrobacter portucalensis", "s__Citrobacter_A amalonaticus", "s__Citrobacter_A farmeri", "s__Citrobacter_A telavivensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter freundii" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter freundii / braakii" & sourmash_species %in% c("s__Citrobacter europaeus", "s__Citrobacter freundii", "s__Citrobacter portucalensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter freundii / braakii" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter freundii complex" & sourmash_species %in% c("s__Citrobacter freundii", "s__Citrobacter portucalensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter freundii complex" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter koseri" & sourmash_species %in% c("Citrobacter koseri", "s__Citrobacter_B koseri"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter koseri" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter koseri-amalonaticus/farmeri" & sourmash_species %in% c("s__Citrobacter_A telavivensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter koseri-amalonaticus/farmeri" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter murliniae" & sourmash_species %in% c("s__Citrobacter murliniae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter murliniae" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter sedlakii" & sourmash_species %in% c("s__Citrobacter_A sedlakii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter sedlakii" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter spp." & sourmash_species %in% c("s__Citrobacter freundii", "Citrobacter freundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter spp." & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter youngae" & sourmash_species %in% c("s__Citrobacter freundii", "s__Citrobacter youngae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Citrobacter youngae" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    # CORYNEBACTERIUM
    mutate(contaminated=ifelse(MALDI_species=="Corynebacterium striatum" & sourmash_species %in% c("s__Corynebacterium striatum"), FALSE, contaminated)) %>%
    # DELFTIA
    mutate(contaminated=ifelse(MALDI_species=="Delftia acidovorans" & sourmash_species %in% c("s__Comamonas tsuruhatensis"), FALSE, contaminated)) %>%
    # ENTEROBACTER
    # Kind of same problem as for  Citrobacter, species level isn't well discriminated
    # E. aerogenes no problem ; E. asburiae (~)
    # E. cloacae includes almost anything from Enterobacter genus
    # Some E. cloacae (n=16), aren't classified by Sourmash but are classified as cronobacter by MLST => don't know what to think about this. Cronobacter should be a separate clade from Enterobacter samples
    # AH Il s'en sort quand même bien avec : E. kobei, E. roggenkampii
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter aerogenes" & sourmash_species %in% c("Klebsiella aerogenes", "s__Klebsiella aerogenes"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter aerogenes" & is.na(sourmash_species) & MLST_species %in% c("kaerogenes"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter asburiae" & sourmash_species %in% c("s__Enterobacter asburiae", "s__Enterobacter asburiae_A", "s__Enterobacter asburiae_B", "s__Enterobacter chengduensis", "s__Enterobacter cloacae_M"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter asburiae" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae" & sourmash_species %in% c("Enterobacter asburiae", "Enterobacter cloacae", "Enterobacter ludwigii", "s__Enterobacter asburiae", "s__Enterobacter asburiae_B", "s__Enterobacter bugandensis", "s__Enterobacter cloacae", "s__Enterobacter hormaechei", "s__Enterobacter hormaechei_A", "s__Enterobacter kobei", "s__Enterobacter ludwigii", "s__Enterobacter quasihormaechei", "s__Enterobacter roggenkampii", "s__Enterobacter sichuanensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae complex" & sourmash_species %in% c("Enterobacter cloacae", "s__Enterobacter asburiae", "s__Enterobacter asburiae_B", "s__Enterobacter bugandensis", "s__Enterobacter chengduensis", "s__Enterobacter cloacae", "s__Enterobacter cloacae_M", "s__Enterobacter cloacae_O", "s__Enterobacter hormaechei_A", "s__Enterobacter kobei", "s__Enterobacter ludwigii", "s__Enterobacter quasihormaechei", "s__Enterobacter roggenkampii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae complex" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae ssp. cloacae" & sourmash_species %in% c("s__Enterobacter cloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae ssp. cloacae" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae ssp. dissolvens" & sourmash_species %in% c("s__Enterobacter hormaechei_A"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter cloacae ssp. dissolvens" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter hormaechei" & sourmash_species %in% c("s__Enterobacter hormaechei_A", "s__Enterobacter quasihormaechei"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter hormaechei" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter hormaechei ssp. hoffmannii" & sourmash_species %in% c("s__Enterobacter hormaechei_A"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter hormaechei ssp. hoffmannii" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter kobei" & sourmash_species %in% c("s__Enterobacter kobei"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter kobei" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter roggenkampii" & sourmash_species %in% c("s__Enterobacter roggenkampii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterobacter roggenkampii" & is.na(sourmash_species) & MLST_species %in% c("ecloacae"), FALSE, contaminated)) %>%
    # ENTEROCOCCUS
    mutate(contaminated=ifelse(MALDI_species=="Enterococcus faecalis" & sourmash_species %in% c("s__Enterococcus faecalis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterococcus faecalis" & is.na(sourmash_species) & MLST_species %in% c("efaecalis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterococcus faecium" & sourmash_species %in% c("Enterococcus faecium", "s__Enterococcus_B faecium"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Enterococcus faecium" & is.na(sourmash_species) & MLST_species %in% c("efaecium"), FALSE, contaminated)) %>%
    # ESCHERICHIA
    mutate(contaminated=ifelse(MALDI_species=="Escherichia coli" & sourmash_species %in% c("Escherichia coli", "s__Escherichia coli", "s__Escherichia marmotae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Escherichia coli" & is.na(sourmash_species) & MLST_species %in% c("ecoli", "ecoli_achtman_4", "escherichia"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Escherichia fergusonii" & sourmash_species %in% c("s__Escherichia fergusonii"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Escherichia fergusonii" & is.na(sourmash_species) & MLST_species %in% c("escherichia"), FALSE, contaminated)) %>%
    # HAEMOPHILUS
    mutate(contaminated=ifelse(MALDI_species=="Haemophilus influenzae" & sourmash_species %in% c("Haemophilus influenzae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Haemophilus influenzae" & is.na(sourmash_species) & MLST_species %in% c("hinfluenzae"), FALSE, contaminated)) %>%
    # HAFNIA
    mutate(contaminated=ifelse(MALDI_species=="Hafnia alvei" & sourmash_species %in% c("Hafnia alvei", "s__Hafnia alvei", "s__Hafnia paralvei", "s__Hafnia proteus"), FALSE, contaminated)) %>%
    # KLEBSIELLA
    mutate(contaminated=ifelse(MALDI_species=="Klebsiella pneumoniae" & sourmash_species %in% c("Klebsiella pneumoniae", "s__Klebsiella africana", "s__Klebsiella pneumoniae", "s__Klebsiella quasipneumoniae", "s__Klebsiella variicola"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Klebsiella pneumoniae" & is.na(sourmash_species) & MLST_species %in% c("klebsiella", "kpneumoniae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Klebsiella oxytoca" & sourmash_species %in% c("Klebsiella michiganensis", "Klebsiella oxytoca", "s__Klebsiella grimontii", "s__Klebsiella michiganensis", "s__Klebsiella oxytoca", "s__Klebsiella planticola", "Klebsiella sp. OBRC7"), FALSE, contaminated)) %>% 
    mutate(contaminated=ifelse(MALDI_species=="Klebsiella oxytoca" & is.na(sourmash_species) & MLST_species=="koxytoca", FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Klebsiella variicola" & sourmash_species=="s__Klebsiella variicola", FALSE, contaminated)) %>%
    # LACTOBACILLUS
    mutate(contaminated=ifelse(MALDI_species=="Lactobacillus casei / paracasei / rhamnosus" & sourmash_species %in% c("s__Lacticaseibacillus rhamnosus"), FALSE, contaminated)) %>%
    # MORAXELLA
    mutate(contaminated=ifelse(MALDI_species=="Moraxella catharralis" & sourmash_species %in% c("s__Moraxella catarrhalis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Moraxella catharralis" & is.na(sourmash_species) & MLST_species %in% c("mcatarrhalis_achtman_6"), FALSE, contaminated)) %>%
    # MORGANELLA
    mutate(contaminated=ifelse(MALDI_species=="Morganella morganii" & sourmash_species %in% c("Morganella morganii", "s__Morganella morganii", "s__Morganella morganii_A"), FALSE, contaminated)) %>%
    # NEISSERIA
    mutate(contaminated=ifelse(MALDI_species=="Neisseria gonorrhoeae" & sourmash_species %in% c("s__Neisseria gonorrhoeae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Neisseria gonorrhoeae" & is.na(sourmash_species) & MLST_species %in% c("neisseria"), FALSE, contaminated)) %>%
    # PASTEURELLA 
    mutate(contaminated=ifelse(MALDI_species=="Pasteurella multocida" & sourmash_species %in% c("s__Pasteurella multocida"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pasteurella multocida" & is.na(sourmash_species) & MLST_species %in% c("pmultocida"), FALSE, contaminated)) %>%
    # PROPIONIBACTERIUM
    mutate(contaminated=ifelse(MALDI_species=="Propionibacterium acnes" & sourmash_species %in% c("s__Cutibacterium acnes"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Propionibacterium acnes" & is.na(sourmash_species) & MLST_species %in% c("pacnes", "pacnes_3"), FALSE, contaminated)) %>%
    # PROTEUS
    mutate(contaminated=ifelse(MALDI_species=="Proteus mirabilis" & sourmash_species %in% c("Proteus mirabilis", "s__Proteus mirabilis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Proteus penneri" & sourmash_species %in% c("Proteus mirabilis", "s__Proteus sp003144375"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Proteus vulgaris" & sourmash_species %in% c("Proteus vulgaris", "s__Proteus sp003144375"), FALSE, contaminated)) %>%
    # PSEUDOCITROBACTER
    # I kinda miss knowledge on this, and unique sample seems to be contaminated
    mutate(contaminated=ifelse(MALDI_species=="Pseudocitrobacter faecalis" & is.na(sourmash_species) & MLST_species %in% c("cfreundii"), FALSE, contaminated)) %>%
    # PSEUDOMONAS
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas aeruginosa" & sourmash_species %in% c("Pseudomonas aeruginosa", "s__Pseudomonas aeruginosa", "s__Pseudomonas_E peradeniyensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas aeruginosa" & is.na(sourmash_species) & MLST_species %in% c("paeruginosa"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas citronellolis" & sourmash_species %in% c("s__Pseudomonas citronellolis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas citronellolis" & is.na(sourmash_species) & MLST_species %in% c("paeruginosa"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas fluorescens" & sourmash_species %in% c("s__Pseudomonas_E fluorescens_BV", "s__Pseudomonas_E lactis", "s__Pseudomonas_E moraviensis_A", "s__Pseudomonas_E neuropathica", "s__Pseudomonas_E paracarnis", "s__Pseudomonas_E sp003033885"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas fluorescens" & is.na(sourmash_species) & MLST_species %in% c("pfluorescens"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas mosselii" & sourmash_species %in% c("s__Pseudomonas_E soli"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas mosselii" & is.na(sourmash_species) & MLST_species %in% c("pputida"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas putida" & sourmash_species %in% c("s__Pseudomonas_E juntendi", "s__Pseudomonas_E ceruminis", "s__Pseudomonas_E fulva", "s__Pseudomonas_E peradeniyensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Pseudomonas putida" & is.na(sourmash_species) & MLST_species %in% c("pputida"), FALSE, contaminated)) %>%
    # RAOULTELLA
    mutate(contaminated=ifelse(MALDI_species=="Raoultella ornithinolytica" & sourmash_species %in% c("s__Klebsiella planticola"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Raoultella ornithinolytica" & is.na(sourmash_species) & MLST_species %in% c("koxytoca"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Raoultella planticola" & sourmash_species %in% c("s__Klebsiella planticola"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Raoultella planticola" & is.na(sourmash_species) & MLST_species %in% c("koxytoca"), FALSE, contaminated)) %>%
    # RHIZOBIUM
    mutate(contaminated=ifelse(MALDI_species=="Rhizobium radiobacter" & sourmash_species %in% c("s__Agrobacterium pusense", "s__Agrobacterium sp900013535", "s__Agrobacterium fabrum"), FALSE, contaminated)) %>%
    # SERRATIA
    # Serratia species are not well discriminated
    mutate(contaminated=ifelse(MALDI_species=="Serratia marcescens" & sourmash_species %in% c("s__Serratia bockelmannii", "s__Serratia marcescens_K", "s__Serratia nevei", "Serratia marcescens", "Serratia sp. HMSC15F11"), FALSE, contaminated)) %>%
    # SHIGELLA
    # Shigella are close to Escherichia genomically also
    mutate(contaminated=ifelse(MALDI_species=="Shigella sonnei" & sourmash_species %in% c("s__Escherichia coli"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Shigella sonnei" & is.na(sourmash_species) & MLST_species %in% c("ecoli_achtman_4"), FALSE, contaminated)) %>%
    # STAPHYLOCOCCUS
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus aureus" & sourmash_species %in% c("s__Staphylococcus aureus", "Staphylococcus aureus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus aureus" & is.na(sourmash_species) & MLST_species %in% c("saureus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus capitis" & sourmash_species %in% c("s__Staphylococcus capitis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus caprae" & sourmash_species %in% c("s__Staphylococcus caprae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus epidermidis" & sourmash_species %in% c("s__Staphylococcus epidermidis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus epidermidis" & is.na(sourmash_species) & MLST_species %in% c("sepidermidis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus haemolyticus" & sourmash_species %in% c("s__Staphylococcus haemolyticus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus haemolyticus" & is.na(sourmash_species) & MLST_species %in% c("shaemolyticus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus hominis" & sourmash_species %in% c("s__Staphylococcus hominis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus hominis" & is.na(sourmash_species) & MLST_species %in% c("shominis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus lugdunensis" & sourmash_species %in% c("s__Staphylococcus lugdunensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus lugdunensis" & is.na(sourmash_species) & MLST_species %in% c("staphlugdunensis"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus pettenkoferi" & sourmash_species %in% c("s__Staphylococcus pettenkoferi"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus saprophyticus" & sourmash_species %in% c("s__Staphylococcus saprophyticus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Staphylococcus simulans" & sourmash_species %in% c("s__Staphylococcus simulans"), FALSE, contaminated)) %>%
    # STENOTROPHOMONAS
    mutate(contaminated=ifelse(MALDI_species=="Stenotrophomonas maltophilia" & sourmash_species %in% c("s__Stenotrophomonas maltophilia", "s__Stenotrophomonas maltophilia_A", "s__Stenotrophomonas maltophilia_AJ", "s__Stenotrophomonas maltophilia_G", "s__Stenotrophomonas maltophilia_O", "Stenotrophomonas maltophilia"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Stenotrophomonas maltophilia" & is.na(sourmash_species) & MLST_species %in% c("smaltophilia"), FALSE, contaminated)) %>%
    # STREPTOCOCCUS
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus anginosus" & sourmash_species %in% c("s__Streptococcus anginosus_C", "s__Streptococcus anginosus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus constellatus" & sourmash_species %in% c("s__Streptococcus constellatus"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus dysgalactiae" & sourmash_species %in% c("s__Streptococcus dysgalactiae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus dysgalactiae" & is.na(sourmash_species) & MLST_species %in% c("sdysgalactiae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus pneumoniae" & sourmash_species %in% c("s__Streptococcus pneumoniae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus pneumoniae" & is.na(sourmash_species) & MLST_species %in% c("spneumoniae"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus pyogenes" & sourmash_species %in% c("s__Streptococcus pyogenes", "Streptococcus pyogenes"), FALSE, contaminated)) %>%
    mutate(contaminated=ifelse(MALDI_species=="Streptococcus pyogenes" & is.na(sourmash_species) & MLST_species %in% c("spyogenes"), FALSE, contaminated)) %>%
    # TURICELLA
    mutate(contaminated=ifelse(MALDI_species=="Turicella otitidis" & sourmash_species %in% c("s__Corynebacterium otitidis"), FALSE, contaminated)) %>%
    
    arrange(genome)
  # Writing to a file
  fwrite(study_df, file = "assembly_summary_illumina_only.tsv", sep = "\t")
  today_date <- Sys.Date()
  # Strains that pass QC
  cat(paste(study_df[NEW_QC_illumina_pass == TRUE]$id, collapse = "\n"), file = paste(today_date,"illumina_only_QC_passed.txt"))
  # Strains that fail Illumina QC
  cat(paste(study_df[NEW_QC_illumina_pass == FALSE]$id, collapse = "\n"), file = paste(today_date,"illumina_only_QC_illumina_failed.txt"))
  
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
     pilon_corrections, pilon_coverage, s, x, nano_stats, n50_stats, line_splitted, line_of_interest,
     line_with_n_reads, N50, N90, nano_mean_read_length, nano_median_read_length, sourmash_species)
}
}

rm(illumina_samples, nano_illumina_samples)
