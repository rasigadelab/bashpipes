#!/bin/bash

rm -r .nextflow
rm -r work
rm .n*
# rm report-*
# rm trace-*

dos2unix --allow-chown src/workflows/workflow.nf
dos2unix --allow-chown src/modules/*.nf
dos2unix --allow-chown src/main.nf
dos2unix --allow-chown src/nextflow.config
dos2unix --allow-chown ./launch_pipeline.sh
