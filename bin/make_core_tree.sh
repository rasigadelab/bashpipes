#!/bin/bash

# Compute fast ML tree from core-genome alignment
# prokaa > panaroo > iqtree

cd sequences

# Fast prokka aligment
launch2.sh 10 sub_prokka samples.txt

# once prokka has run, gather GFF files for panaroo 

cd ..
mkdir -p panaroo/gff

cp ./sequences/*/prokka/*.gff panaroo/gff 
cd panaroo
panaroo -i ./gff/*.gff -o . --clean-mode strict -a core --core_threshold 1 -t 24 

cd ..
mkdir iqtree  
cp panaroo/core_gene_alignment.aln iqtree 
cd iqtree 

iqtree -s core_gene_alignment.aln -m HKY