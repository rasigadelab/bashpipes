#!/bin/bash
# PROKKA Annotation

SAMPLE=$1
SAMPLE_POLISHED="$SAMPLE""_polished_realigned"

cd $SAMPLE

mkdir -p mlst
mkdir -p sourmash
mkdir -p amrfinder
mkdir -p prokka
mkdir -p mob_recon

# MLST
# Determine sequence type (ST) of the isolate
cd mlst &&
mlst ../"$SAMPLE_POLISHED".fasta > mlst.tsv 2>> mlst.log
cd ..

# TAXON VERIFICATION
# Use Rscript to extract taxon from sourmash output. This writes files genus.txt and species.txt
# TODO formalize database location

# ./SOURMASH directory
cd sourmash &&
sourmash sketch dna -p scaled=1000,k=21 ../"$SAMPLE_POLISHED".fasta &> sourmash.log &&
sourmash lca classify --query "$SAMPLE_POLISHED".fasta.sig --db ~/common_db/sourmash-genbank-k21.json.gz > sourmash.csv 2>> sourmash.log &&
sourmash_extract.R

# STRAIN WORKING DIRECTORY
cd ..


# Extract genus and species names if available
GENUS=""

if [ -f genus.txt ]; then
    GENUS=$(<genus.txt)
    # echo $GENUS
fi

SPECIES=""

if [ -f species.txt ]; then
    SPECIES=$(<species.txt)
    # echo $SPECIES
fi


# AMRFINDER

# ORGANISM-DEPENDENT AMRFINDER
# is picky about organism so add flag only if valid organism

# amrfinder --list_organisms
# Acinetobacter_baumannii, Burkholderia_cepacia, Burkholderia_pseudomallei, Campylobacter, Clostridioides_difficile, Enterococcus_faecalis, Enterococcus_faecium, Escherichia, Klebsiella, Neisseria, Pseudomonas_aeruginosa, Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius, Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae

ORGANISM_FLAG=""

case "$GENUS" in
    "Klebsiella" ) ORGANISM_FLAG="--organism Klebsiella";;
    "Escherichia" ) ORGANISM_FLAG="--organism Escherichia";;
    *) ORGANISM_FLAG=""
esac

# echo $ORGANISM_FLAG

amrfinder --plus --nucleotide "$SAMPLE_POLISHED".fasta --threads 32 $ORGANISM_FLAG > amrfinder/amrfinder.tsv 2> amrfinder/amrfinder.log

# ANNOTATION

# conda activate issue, see https://github.com/conda/conda/issues/7980

# STRAIN DIRECTORY

# cd .. &&
source ~/miniconda3/etc/profile.d/conda.sh &&
conda activate prokka &&
prokka "$SAMPLE_POLISHED".fasta --force –-addgenes –-cpus 0 \
    --outdir prokka --compliant --prefix $SAMPLE \
    –-usegenus –-genus $GENUS –-species $SPECIES \
    &> prokka/prokka.log
conda deactivate

## MGE ANALYSIS
# cd .. &&
mob_recon -n 32 --force --infile "$SAMPLE_POLISHED".fasta --outdir mob_recon \
    1> mob_recon/mob_recon.log 2> mob_recon/mob_recon.err

echo "DONE" > annotate_done.flag