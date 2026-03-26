#!/bin/bash
#
# Title: dispatch.sh
# Description: Creates a sample directory and copies sample reads inside it.
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

while read SAMPLE; do
    mkdir $SAMPLE &&
    mv "$SAMPLE"*.fastq.gz ./"$SAMPLE"
done < ./samples.txt