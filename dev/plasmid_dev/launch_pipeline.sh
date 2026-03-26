#!/bin/bash
#
# Title: launch_pipeline.sh
# Description: Runs plasmid comparison pipeline.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

nextflow -C src/nextflow.config run src/main.nf -params-file src/params_plasmid_compa.json -with-trace -with-report --prefix run_plasmid -profile standard

