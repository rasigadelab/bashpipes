#!/bin/bash
#
# Title: launch_pipeline.sh
# Description: Runs phylogenetic pipelines, first MASH, then Variant Calling and finally Phylogeny.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_mash.json -with-trace -with-report --prefix run_1 -profile standard

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_variant_calling.json -with-trace -with-report --prefix run_1 -profile standard

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_phylogeny.json -with-trace -with-report --prefix run_1 -profile standard
