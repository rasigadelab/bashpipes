#!/bin/bash

nextflow -C src/nextflow.config run src/main.nf -params-file src/params_plasmid_compa.json -with-trace -with-report --prefix run_plasmid -profile standard

