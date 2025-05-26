#!/bin/bash

nextflow -C src/nextflow.config run src/main.nf -params-file src/params_mash.json -with-trace -with-report --prefix run_1 -profile standard

# nextflow -C src/nextflow.config run src/main.nf -params-file src/params_variant_calling.json -with-trace -with-report --prefix run_1 -profile standard

# nextflow -C src/nextflow.config run src/main.nf -params-file src/params_phylogeny.json -with-trace -with-report --prefix run_1 -profile standard
