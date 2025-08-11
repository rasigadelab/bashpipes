#!/bin/bash

NXF_OFFLINE=true NXF_DISABLE_AUTO_UPDATE=true nextflow -C src/nextflow.config run src/main.nf -params-file src/params_nano_illumina.json -with-trace -with-report --prefix run_nano_illumina -profile standard

