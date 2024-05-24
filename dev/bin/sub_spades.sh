#!/bin/bash

SAMPLE=$1

cd $SAMPLE
mkdir spades

spades.py -1 "$SAMPLE"_R1.fastq.gz -2 "$SAMPLE"_R2.fastq.gz --isolate -o spades \
    1> spades/spades_1.log 2> spades/spades.err

cp spades/contigs.fasta "$SAMPLE""_polished.fasta"