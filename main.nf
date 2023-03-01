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

process preprocessFASTA {

    container "${params.container__biopython}"
    label "io_limited"
    
    input:
    tuple val(sample_name), file(fasta)

    output:
    tuple val(sample_name), file("${sample_name}.fasta")

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
    tuple val(sample_name), file(fasta)

    output:
    path "${sample_name}*"

"""
#!/bin/bash

set -euxo pipefail

echo Running Prokka

prokka \
    --outdir "OUTPUT/" \
    --prefix "${sample_name}" \
    --cpus ${task.cpus} \
    "${fasta}"

echo Compressing outputs

mv OUTPUT/* ./
rmdir OUTPUT
gzip "${sample_name}"*

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
    tuple val(sample_name), file(fasta)

    output:
    path "*"

"""
#!/bin/bash

set -euxo pipefail

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

echo Running gapseq

gapseq doall "${fasta}"

echo Done

"""

}

process collect_gapseq {

    container "${params.container__widgets}"
    label "io_limited"
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
    errorStrategy 'retry'
    maxRetries 2

    input:
    path "*"

    output:
    path "*.html"

"""
gapseq-widget.py
"""

}


workflow {

    // Print the help message
    if (params.help){
        helpMessage();
        exit 0
    }

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
        .from(file(params.sample_sheet, checkIfExists: true))
        .splitCsv(header:true)
        .map { it -> [it["name"], file(it["fasta"])]}
        .set { sample_sheet_ch }

    preprocessFASTA(
        sample_sheet_ch
    )

    run_prokka(preprocessFASTA.out)
    run_gapseq(preprocessFASTA.out)

    collect_gapseq(
        run_gapseq.out.flatten().toSortedList()
    )

}
