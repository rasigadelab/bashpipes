#!/bin/bash

# Gather mob_recon output files into some new directory, mimicking file structure

NEWDIR=$1/sequences

KEY=$2

# Find all matching files
FNAMES=$(ls ./*/mob_recon/$KEY.fasta) 

# Find corresponding isolate names (these should be directory names)
INAMES=$(echo $FNAMES | grep -Po "(?<=./)([^/^\s]+)(?=/mob_recon)") 

mkdir -p $NEWDIR

for INAME in $INAMES; do mkdir $NEWDIR/$INAME; done 

for INAME in $INAMES; do cp ./$INAME/mob_recon/$KEY.fasta $NEWDIR/$INAME; done

# Save default fofn, replacing spaces with newlines as required by launch2 format
echo $INAMES | tr ' ' '\n' > $NEWDIR/samples.txt

#cat samples.txt | tr ' ' '\n' > samples.txt