#!/bin/bash

rm -r .nextflow
rm -r work
rm -r ./workflows/.nextflow
rm -r ./workflows/work
rm ./workflows/.n*
rm .n* 

dos2unix ./workflows/workflow.nf
dos2unix ./modules/*.nf
dos2unix ./main.nf
dos2unix ./nextflow.config
dos2unix ./launch_test.sh