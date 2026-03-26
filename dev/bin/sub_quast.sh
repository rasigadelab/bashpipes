#!/bin/bash
#
# Title: sub_quast.sh
# Description: Launches QUAST, QC reporting on polished genome assembly.
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

SAMPLE=$1

cd $SAMPLE
mkdir quast

quast.py -o quast "$SAMPLE""_polished.fasta" \
    1> quast/quast.log \
    2> quast/quast.err