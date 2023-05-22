#!/bin/bash

 nextflow -C nextflow.config run main.nf -params-file params_phylogeny.json -with-trace -with-report --prefix run_1 -profile developer -stub-run
   