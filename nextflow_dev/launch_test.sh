#!/bin/bash

 nextflow -C nextflow.config run main.nf -params-file params_illumina_only.json -stub-run -with-trace --prefix run_1 -profile developer
   