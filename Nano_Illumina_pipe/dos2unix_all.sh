#!/bin/bash
#
# Title: dos2unix_all.sh
# Description: Cleans nextflow logs (.nextflow file and work directory) and convert to unix format scripts to execute.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

rm -r .nextflow
rm -r work
rm .n*
# rm report-*
# rm trace-*

dos2unix --allow-chown src/workflows/workflow.nf
dos2unix --allow-chown src/modules/*.nf
dos2unix --allow-chown src/main.nf
dos2unix --allow-chown src/nextflow.config
dos2unix --allow-chown ./launch_pipeline.sh
