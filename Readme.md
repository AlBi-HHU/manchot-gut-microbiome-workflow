Snakemake workflow for analysis of metagenomic microbiome data in either fastq or fast5 format

# Setup

Requires Snakemake and Conda available on the system

# Execution

## Required Input

- Metagenomic Samples as fast5/fastq or barcoded fast5
- Human Reference Genome (for filtering purposes)
- Suitable metagenomics .fa database 

Additional input such as illumina sequencing based reads as well as additional databases are required for optional components of the workflow. Refer to the config.yaml for details.

## Configuration

Move (or symlink) files and folders to data/input such that

- data/input/databases should contain the databases 
- data/input/reads contains fastq files
- data/input/signals contains fast5 files

Copy config.example.yaml, samples.example.yaml and signals.example.yaml.

Remove the example infix from the filenames, adjust file contents if needed

Individual components of the workflow can be toggled via the config.yaml file

# Notebook

A jupyter notebook with downstream analysis and scripts to create visualizations used in the paper is included in the notebook folder.