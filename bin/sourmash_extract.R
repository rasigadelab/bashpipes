#!/usr/bin/env Rscript

# Extract genus and species from sourmash output
tb <- read.csv("sourmash.csv")

genus <- tb$genus
if(!is.na(genus)) cat(genus, file = "genus.txt")

species <- tb$species
if(!is.na(species)) {
    species <- strsplit(species, " ")[[1]][2]
    cat(species, file = "species.txt")
}