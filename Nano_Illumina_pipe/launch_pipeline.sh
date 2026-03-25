#!/bin/bash
#
# Title: launch_pipeline.sh
# Description: Runs Nano Illumina assembly/annotation pipeline.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_nano_illumina.json -with-trace -with-report --prefix run_nano_illumina -profile standard

