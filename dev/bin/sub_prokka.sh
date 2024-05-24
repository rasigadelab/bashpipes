#!/bin/bash

# Run fast prokka in parallel on subdirectories

SAMPLE=$1

cd $SAMPLE
mkdir prokka

source ~/miniconda3/etc/profile.d/conda.sh &&
conda activate prokka &&
prokka *.fasta --force –-addgenes –-cpus 0 \
    --outdir prokka --compliant --prefix $SAMPLE \
    &> prokka/prokka.log
conda deactivate