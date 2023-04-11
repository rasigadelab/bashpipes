#!/bin/bash

 nextflow -C nextflow.config run main.nf -params-file params_illumina_only.json -with-trace -with-report --prefix run_1 -profile standard
   