#!/bin/bash

rm -r .nextflow
rm -r work
rm .n*
#rm trace-*
#rm report-*

dos2unix ./src/workflows/workflow.nf
dos2unix ./src/modules/*.nf
dos2unix ./src/main.nf
dos2unix ./src/nextflow.config
dos2unix ./launch_pipeline.sh
