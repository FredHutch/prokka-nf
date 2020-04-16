#!/bin/bash

set -e

NXF_VER=20.01.0 nextflow \
    -c ~/nextflow.config \
    run \
    main.nf \
    --sample_sheet test/sample_sheet.csv \
    --output_folder test/output \
    -resume \
    -with-trace
