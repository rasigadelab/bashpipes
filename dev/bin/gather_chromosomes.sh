#!/bin/bash
# Gather all chromosome fasta after mob_recon run

mkdir -p molecules

for SAMPLE in $(<samples.txt);
 do
    echo $SAMPLE
    FASTA=$(ls "$SAMPLE"/mob_recon/*.fasta)
    # FASTA<ls "$SAMPLE"/mob_recon/chromosome.fasta
    for CHROM in $FASTA;
    do
        # Prepend sample name to copied assembly
        TARGET=${CHROM/"$SAMPLE""/mob_recon/"/}
        TARGET="molecules/""$SAMPLE"_"$TARGET"
        cp $CHROM $TARGET

        # Prepend sample name to fasta headers
        sed -i "s/^>/>${SAMPLE}_/g" $TARGET

        # echo $TARGET
        # TARGET=${CHROM/chromosome/"$SAMPLE"_chromosome}
        # echo $TARGET
        # TARGET=${TARGET/'"$SAMPLE"/mob_recon'/'chromosomes/'}
    done
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