#!/bin/bash
#
# Title: mob_recon_gather.sh
# Description: Launches Prokka annotation, Panaroo (pan-genome analysis) and computes fast Maximum-Likelihood tree (IQtree) from core-genome alignment
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

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