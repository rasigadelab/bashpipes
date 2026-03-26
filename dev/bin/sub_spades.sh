#!/bin/bash
#
# Title: sub_spades.sh
# Description: Launches SPAdes assembler on R1/R2 Illumina reads.
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

SAMPLE=$1

cd $SAMPLE
mkdir spades

spades.py -1 "$SAMPLE"_R1.fastq.gz -2 "$SAMPLE"_R2.fastq.gz --isolate -o spades \
    1> spades/spades_1.log 2> spades/spades.err

cp spades/contigs.fasta "$SAMPLE""_polished.fasta"