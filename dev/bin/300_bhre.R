rm(list = ls())

library(data.table)
library(readxl)
library(stringr)
library(openxlsx)

# HELP TO CHECK QC FOR BIG BATCH OF DATA
metadata_df <- data.table(read_excel("metadata.xlsx", sheet = 1))
sample_list <- data.table(read.table("all_epitrack_to_assemble.txt", header = FALSE, sep = ";"))
colnames(sample_list) <- c("SAMPLE_ID", "CLUSTER_ID")

tmp <- metadata_df[SAMPLE_ID %in% sample_list$SAMPLE_ID,]
tmp <- tmp[,c("SAMPLE_ID", "ordre chronologique", "SPECIES", "SEQ_ILLUMINA_DATE")]
tmp[,"Nb_Illumina_seq"] <- ifelse(is.na(tmp$SEQ_ILLUMINA_DATE), 1, str_count(tmp$SEQ_ILLUMINA_DATE, "/")+1)
colnames(tmp) <- c("SAMPLE_ID", "ORDRE", "MALDI_species", "SEQ_ILLUMINA_DATE", "NB_ILLUMINA_SEQ")
# GO CHECK QC file for QC_PASS + Sourmash species + MLST species
qc_df <- data.table(read.table("./genomes/assembly_summary_illumina_only.tsv", header = TRUE, sep = "\t"))
qc_df <- qc_df[,c("id", "species", "sourmash_species", "QC_illumina_pass")]
colnames(qc_df) <- c("id", "MLST_species", "SOURMASH_species", "QC_illumina_pass")

# Merge 
tmp <- merge(tmp, qc_df, by.x = "SAMPLE_ID", by.y = "id")
rm(metadata_df, qc_df, sample_list)

# Contamination track
tmp2 <- tmp
tmp2$SOURMASH_species <- gsub("s__", "", tmp$SOURMASH_species)
tmp2$SOURMASH_species <- str_to_lower(paste0(substr(tmp2$SOURMASH_species, 1, 1), str_extract(tmp2$SOURMASH_species, "[^ ]+$")))
tmp2$MALDI_species <- str_to_lower(paste0(substr(tmp$MALDI_species, 1, 1), str_extract(tmp$MALDI_species, "[^ ]+$")))
citrobacter <- c("cbraakii", "cfarmeri", "cfreundii", "ckoseri", "csedlakii", "cspp", "cyoungae", "ccomplex")
enterobacter <- c("ecloacae", "ecomplex", "ehormaechei", "espp.", "easburiae", "eroggenkampii",
                  "ekobei", "equasihormaechei", "ehormaechei_a", "emarmotae")
ecoli <- c("ecoli_achtman_4")
kpneumoniae <- c("klebsiella", "kpneumoniae", "kquasipneumoniae", "kvariicola", "kafricana")
kaerogenes <- c("eaerogenes", "kaerogenes")
koxytoca <- c("koxytoca", "kmichiganensis", "kgrimontii", "kplanticola")

redef_species <- function(list_of_annot, new_species){
  tmp2$SOURMASH_species <- ifelse(tmp2$SOURMASH_species %in% list_of_annot, new_species, tmp2$SOURMASH_species)
  tmp2$MALDI_species <- ifelse(tmp2$MALDI_species %in% list_of_annot, new_species, tmp2$MALDI_species)
  tmp2$MLST_species <- ifelse(tmp2$MLST_species %in% list_of_annot, new_species, tmp2$MLST_species)
  return(tmp2)
}
tmp2 <- redef_species(enterobacter, "enterobacter")
tmp2 <- redef_species(citrobacter, "citrobacter")
tmp2 <- redef_species(kpneumoniae, "kpneumoniae")
tmp2 <- redef_species(ecoli, "ecoli")
tmp2 <- redef_species(kaerogenes, "kaerogenes")
tmp2 <- redef_species(koxytoca, "koxytoca")

# maybe add here later, the case with missing info
tmp2$CONTA_check <- ifelse(mapply(grepl,tmp2$MALDI_species, tmp2$SOURMASH_species), "OK",
                           ifelse(tmp2$MALDI_species==tmp2$MLST_species, "OK", "NOK"))
rm(citrobacter, ecoli, enterobacter, kaerogenes, koxytoca, kpneumoniae, redef_species)


# Adding final qc+conta check --> what verifications are needed
recurrent_species <- c("ecoli", "enterobacter", "citrobacter", "paeruginosa", "kpneumoniae", "koxytoca")
final_verif <- c()
# * if QC_PASS is FALSE & (Sourmash species or MLST species is coherent with MALDI species) => 
for(sample in tmp2$SAMPLE_ID){
  current <- tmp2[tmp2$SAMPLE_ID==sample,]
  if(current$QC_illumina_pass==FALSE & current$CONTA_check=="OK"){
    #      if species is ecoli => check manually
    if(current$MALDI_species=="ecoli"){
      current$VERIF <- "MANUAL CHECK REQUIRED"
    } 
    else if(!(current$MALDI_species %in% recurrent_species)){
      current$VERIF <- "MANUAL CHECK REQUIRED"
    }
    else if(current$NB_ILLUMINA_SEQ >=2){
      #      if sample already sequenced twice => sequencing finished
      current$VERIF <- "SEQUENCING FINISHED"
    }
    else{
      #      else => 2nd sequencing required
      current$VERIF <- "2nd SEQUENCING REQUIRED"
    }
  }
  # * if QC_PASS is FALSE & (Sourmash species or MLST species is NOT coherent with MALDI species) =>
  if(current$CONTA_check=="NOK"){
    if(current$NB_ILLUMINA_SEQ >=2){
      #      if sample already sequenced twice => sequencing finished + SAMPLE is NON DISPO
      current$VERIF <- "SAMPLE_NON_DISPO"
    }
    else{
      #      else => MALDI required, 2nd sequencing required if MALDI invalid contamination
      current$VERIF <- "MALDI_REQUIRED AND 2ND SEQUENCING"
    }
  }
  if(current$QC_illumina_pass==TRUE & current$CONTA_check=="OK"){
    current$VERIF <- "ALL GOOD"
  }
  final_verif <- c(final_verif, current$VERIF)
}
tmp3 <- cbind(tmp2, final_verif)
rm(final_verif, recurrent_species, current, sample)

# Output
qc_df <- data.table(read.table("./genomes/assembly_summary_illumina_only.tsv", header = TRUE, sep = "\t"))
qc_df <- qc_df[qc_df$id %in% tmp3$SAMPLE_ID,]
write.xlsx(list(tmp3, qc_df), file = paste0("bhre_checking.xlsx"), rowNames = FALSE, 
           sheetName = c("Checking", "QC_metrics"))

