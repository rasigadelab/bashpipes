#!/bin/bash
#
# Title: launch2.sh
# Description: Launches parallel tasks for assembling (Nano-Illumina).
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#

while read SAMPLE; do
  flyepilon.sh $SAMPLE &
done < ./$1