#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/prokka-nf <args>

    Required Arguments:
        --sample_sheet        # Comma-separated table with `fasta` and `name` columns
        --output_folder       # Folder to place output files

    """.stripIndent()
}

// Print the help message
params.help = false
if (params.help){
    helpMessage();
    exit 0
}

// Default values for flags
params.sample_sheet = false
params.output_folder = false

if (!params.sample_sheet){
    log.info"""Please specify --sample_sheet."""
    exit 0
}
if (!params.output_folder){
    log.info"""Please specify --output_folder."""
    exit 0
}

// Set up a channel with the fasta and yaml paths

Channel
    .from(file(params.sample_sheet))
    .splitCsv(header:true)
    .map { it -> [it["name"], file(it["fasta"])]}
    .set { sample_sheet_ch }

process preprocessFASTA {

    container "${params.container__biopython}"
    label "io_limited"
    
    input:
    tuple sample_name, file(fasta) from sample_sheet_ch

    output:
    tuple val(sample_name), file("${sample_name}.fasta") into for_prokka, for_gapseq

    """
#!/usr/bin/env python3

# Following criteria from https://github.com/ncbi/pgap/wiki/Input-Files

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import re

# Sanitize and write out
def preprocess_fasta(genome, handle):

    seen_headers = set([])
    
    for header, seq in genome.items():
        
        # Make sure the sequence is >= 199 nucleotides
        if len(seq) < 199:
            continue

        # Trim the header to 37 characters
        if len(header) > 37:
            header = header[:37]

        # Only include letters, digits, hyphens (-), underscores (_), periods (.), colons (:), asterisks (*), and number signs (#)
        header = re.sub('[^0-9a-zA-Z-.*#\$_:]', '_', header)

        # All headers are unique
        assert header not in seen_headers, "Duplicate header: %s (note truncation to first 37 characters)" % header
        seen_headers.add(header)

        # Make sure there are no N's at the beginning or end
        assert seq[0] != "#"
        assert seq[-1] != "#"

        handle.write(">%s\\n%s\\n" % (header, seq))

# Read in all of the genome
if "${fasta}".endswith(".gz"):
    genome = dict([
        (header, seq)
        for header, seq in SimpleFastaParser(gzip.open("${fasta}", "rt"))
    ])
else:
    genome = dict([
        (header, seq)
        for header, seq in SimpleFastaParser(open("${fasta}", "r"))
    ])

with open(
    "${sample_name}.fasta", 
    "w"
) as handle:
    preprocess_fasta(genome, handle)



    """

}

process run_prokka {

    container "${params.container__prokka}"
    label "io_limited"
    publishDir "${params.output_folder}/${sample_name}/prokka/", mode: 'copy', overwrite: true
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_name), file(fasta) from for_prokka

    output:
    path "${sample_name}*"

"""
#!/bin/bash

set -euxo pipefail

echo Running Prokka

prokka \
    --outdir "./" \
    --prefix "${sample_name}" \
    --cpus ${task.cpus} \
    "${fasta}"

echo Compressing outputs

gzip "${sample_name}*"

echo Done

"""

}

process run_gapseq {

    container "${params.container__gapseq}"
    label "io_limited"
    publishDir "${params.output_folder}/${sample_name}/gapseq/", mode: 'copy', overwrite: true
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_name), file(fasta) from for_gapseq

    output:
    path "*"

"""
#!/bin/bash

set -euxo pipefail

echo Running gapseq

gapseq doall "${fasta}"

echo Done

"""

}