#!/bin/bash
#
# Title: to_run_if_snippy_find_0_snps.sh
# Description: Commands to run if Snippy process ends with an error because of 0 SNP detected between samples (see tail of snippy.err).
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

# Inside work dir, snippy section
singularity exec --bind $(pwd) /data/containers/snippy.sif snp-dists core.full.aln > snp_matrix.tsv
singularity exec --bind $(pwd) /data/containers/snippy.sif sed -n '/>Reference/q;p' core.full.aln > core_without_ref.aln
singularity exec --bind $(pwd) /data/containers/snippy.sif snippy-clean_full_aln core_without_ref.aln > clean.full.aln
singularity exec --bind $(pwd) /data/containers/snippy.sif snippy-clean_full_aln core.full.aln > clean.full.withref.aln

# Copy results !
cd ..
cp -r snippy /data/Phylogeny_results/phylogeny/ecloacae/minicluster_1

