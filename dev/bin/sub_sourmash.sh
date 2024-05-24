#!/bin/bash

SAMPLE=$1

cd $SAMPLE

mkdir -p sourmash

cd sourmash &&
sourmash sketch dna -p scaled=1000,k=21 ../"$SAMPLE""_polished.fasta" &> sourmash.log &&
sourmash lca classify --query "$SAMPLE""_polished.fasta.sig" --db ~/common_db/sourmash-genbank-k21.json.gz > sourmash.csv 2>> sourmash.log &&
sourmash_extract.R