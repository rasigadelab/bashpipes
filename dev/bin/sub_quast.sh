#!/bin/bash

SAMPLE=$1

cd $SAMPLE
mkdir quast

quast.py -o quast "$SAMPLE""_polished.fasta" \
    1> quast/quast.log \
    2> quast/quast.err