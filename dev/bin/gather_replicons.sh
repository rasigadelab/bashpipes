#!/bin/bash
# Gather all chromosome fasta after mob_recon run
# Specify new directory, mob_recon replicon name (eg chromosome, plasmid_AA002 etc), sample file
# gather_chromosomes_02.sh NEWDIR REPLICON SAMPLES.TXT

NEWDIR=$1/sequences
REPLICON=$2
SAMPLES=$3

mkdir -p $NEWDIR

for SAMPLE in $(<$SAMPLES);
 do
    # echo $SAMPLE
    FASTA=$(ls "$SAMPLE"/mob_recon/"$REPLICON".fasta)

    # Prepend sample name to copied file name
    TARGET="$SAMPLE"_"$REPLICON".fasta

    cp $FASTA $TARGET

    # Prepend sample name to all fasta headers
    sed -i "s/^>/>${SAMPLE}_/g" $TARGET    
 done


# for i in *.faa; do 
#     sed -i "s/^>/>${i}_/g" *.faa
# done

# Example plasmid phylogeny: might be automated by specifying a plasmid cluster ID, here AA002
# mkdir AA002
# cp *AA002.fasta ./AA002
# cd AA002
# cat *.fasta > AA002.fasta
# mafft --thread 32 AA002.fasta > AA002_aligned.fasta
# iqtree2 -s AA002_aligned.fasta

# Alignment is total garbage, just a test. MAFFT very slow in this case (minutes) but IQtree very fast (seconds including model finding).