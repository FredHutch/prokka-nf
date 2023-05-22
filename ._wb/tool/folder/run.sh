#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow from ${PWD}"
echo

# Build the samplesheet
echo "Building the samplesheet"
echo fasta,name > sample_sheet.csv
for fp in ${FOLDER%/}/*${SUFFIX}; do
    n="${fp##*/}"
    n="${n%$SUFFIX}"
    echo "$fp,$n" >> sample_sheet.csv
done


# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}" \
    --output_folder "${PWD}" \
    --sample_sheet sample_sheet.csv \
    -params-file ._wb/tool/params.json \
    -resume

# If temporary files were not placed in a separate location
if [ -d work ]; then
    # Delete the temporary files created during execution
    echo Removing temporary files
    rm -r work
fi

echo
date
echo Done
