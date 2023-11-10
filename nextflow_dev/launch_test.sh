#!/bin/bash

nextflow -C nextflow.config run main.nf -params-file params_nano_illumina.json -with-trace -with-report --prefix run_1 -profile standard

