# Nextflow tool running Prokka

The goal of this tool is to run the Prokka genome annotation tool
on arbitrary sets of genomes using the Nextflow system for workflow management.

## Input Data

In order to run this tool, you will need to provide:

  1) Genome sequences in FASTA format
  2) Sample sheet linking each genome FASTA with its name


### Sample sheet

The sample sheet linking each genome FASTA with its name
must be formatted as a CSV with two named columns, `fasta` and `name`.

For example:

```
fasta,name
local_folder/genome1.fasta,genome1
local_folder/genome2.fasta,genome2
```

The path to the FASTA can be a relative path, an absolute path, or even
the URL of a file which can be accessed via FTP, HTTP, or any other common file protocol.

## Setting up Nextflow

In order to run the tool, you first need to install and configure Nextflow. Take a look
at the [Nextflow documentation](http://nextflow.io/) for help with this. Some instructions
and guidance on this process for Fred Hutch investigators can be found 
[here](https://sciwiki.fredhutch.org/compdemos/nextflow/). 

## Running Prokka

Once you have the genome FASTAs and sample sheet CSV, the last thing you
need to decide is where the output files will be placed. With that in hand you can run
the Prokka pipeline as follows:

```
nextflow \
    run \
    FredHutch/prokka-nf \
    --sample_sheet PATH_TO_SAMPLE_SHEET.CSV \
    --output_folder PATH_TO_OUTPUT_FOLDER/
```

Some elaborations on this command may be provided by, e.g.:

- Specifying a specific version of Nextflow to use (`NXF_VER=20.01.0 nextflow run ...`)
- Picking up from a halted run (`-resume`)
- Specifying a work directory (`-work-dir`)
- Etc.
