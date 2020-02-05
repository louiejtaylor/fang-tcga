# fang-tcga
Snakemake workflow to map non-human TCGA reads to target genomes and summarize alignments

Note: this project is under active development and is not currently operational.

## Install

    git clone https://github.com/louiejtaylor/fang-tcga

Edit the config file to point to your desired output dir and the list of files to 
process (and/or download). Making a new/copy of your config file for each experiment is 
recommended for reproducibility.

## Requirements

Requires Snakemake and Conda to run. All other dependencies are handled at runtime by Conda.

## Run

    snakemake --configfile config.yml --use-conda

If you'd like only certain outputs, specify a target rule like so:

    snakemake --configfile config.yml --use-conda all_preprocess

Available target rules include: `all_download`, `all_preprocess`, `all_align`, `all_summarize`
