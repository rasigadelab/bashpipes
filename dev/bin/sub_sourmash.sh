#!/bin/bash
#
# Title: sub_sourmash.sh
# Description: Launches Sourmash (taxonomic annotation) on polished genome assembly.
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

SAMPLE=$1

cd $SAMPLE

mkdir -p sourmash

cd sourmash &&
sourmash sketch dna -p scaled=1000,k=21 ../"$SAMPLE""_polished.fasta" &> sourmash.log &&
sourmash lca classify --query "$SAMPLE""_polished.fasta.sig" --db ~/common_db/sourmash-genbank-k21.json.gz > sourmash.csv 2>> sourmash.log &&
sourmash_extract.R