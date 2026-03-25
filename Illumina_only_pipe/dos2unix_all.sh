#!/bin/bash
#
# Title: dos2unix_all.sh
# Description: Cleans nextflow logs (.nextflow file and work directory) and convert to unix format scripts to execute.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

rm -rf .nextflow
rm -rf work
rm .n*
# rm report-*
# rm trace-*

dos2unix src/workflows/workflow.nf
dos2unix src/modules/*.nf
dos2unix src/main.nf
dos2unix src/nextflow.config
dos2unix launch_pipeline.sh
