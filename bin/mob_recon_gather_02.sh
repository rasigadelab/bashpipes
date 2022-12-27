#!/bin/bash

# Gather mob_recon output files into some new directory, mimicking file structure

NEWDIR=$1/sequences

REPLICON=$2

SAMPLES=$3

mkdir -p $NEWDIR

for SAMPLE in $(<$SAMPLES)
    do
        # echo $SAMPLE
        FNAME=./"$SAMPLE"/mob_recon/"$REPLICON".fasta
        # echo $FNAME
    
        mkdir -p "$NEWDIR"/"$SAMPLE"

        # Prepend sample name to file name
        TARGET="$SAMPLE"_"$REPLICON".fasta

        cp $FNAME "$NEWDIR"/"$SAMPLE"/"$TARGET"
    done

# Save default fofn, replacing spaces with newlines as required by launch2 format
# echo $SAMPLES | tr ' ' '\n' > $NEWDIR/samples.txt
cp $SAMPLES $NEWDIR/samples.txt

exit



# Find all matching files
FNAMES=$(ls ./*/mob_recon/$REPLICON.fasta) 

# Find corresponding isolate names (these should be directory names)
INAMES=$(echo $FNAMES | grep -Po "(?<=./)([^/^\s]+)(?=/mob_recon)") 

mkdir -p $NEWDIR

for INAME in $INAMES; do mkdir $NEWDIR/$INAME; done 

for INAME in $INAMES; do cp ./$INAME/mob_recon/$REPLICON.fasta $NEWDIR/$INAME; done

# Save default fofn, replacing spaces with newlines as required by launch2 format
echo $INAMES | tr ' ' '\n' > $NEWDIR/samples.txt

#cat samples.txt | tr ' ' '\n' > samples.txt