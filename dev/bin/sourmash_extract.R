#!/usr/bin/env Rscript

# Title: sourmash_extract.R
# Description: Extracts genus and species from Sourmash results and outputs it into genus.txt and species.txt.
# Author: Jean-Philippe Rasigade
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Jean-Philippe Rasigade

# Extract genus and species from sourmash output
tb <- read.csv("sourmash.csv")

genus <- tb$genus
if(!is.na(genus)) cat(genus, file = "genus.txt")

species <- tb$species
if(!is.na(species)) {
    species <- strsplit(species, " ")[[1]][2]
    cat(species, file = "species.txt")
}