#!/bin/bash

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_mash.json -with-trace -with-report --prefix run_1 -profile standard

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_variant_calling.json -with-trace -with-report --prefix run_1 -profile standard

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_phylogeny.json -with-trace -with-report --prefix run_1 -profile standard
