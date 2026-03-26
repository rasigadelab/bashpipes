#!/bin/bash
#
# Title: sub_prokka.sh
# Description: Launches Prokka, global genes predictions on polished genome assembly.
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

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