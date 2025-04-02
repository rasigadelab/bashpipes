#!/bin/bash

rm -rf .nextflow
rm -rf work
rm .n*
# rm report-*
# rm trace-*

dos2unix src/workflows/workflow.nf
dos2unix src/modules/*.nf
dos2unix src/main.nf
dos2unix src/nextflow.config
dos2unix launch_pipeline.sh
