# fang-tcga
Snakemake workflow to map TCGA reads to target genomes and summarize alignments

Warning: this project is under active development and may break in unexpected ways.

## Install

    git clone https://github.com/louiejtaylor/fang-tcga

Edit the config file to point to your desired output dir and the list of files to 
process (and/or download). Making a new/copy of your config file for each experiment is 
recommended for reproducibility.

## Requirements

Requires Snakemake and Conda to run. All other dependencies are handled at runtime by Conda.

## Run

**Note: if the data you are trying to access are controlled (e.g., not public), you will need an access token from the Genomic Data Commons.** There are detailed instructions for obtaining one [here](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/#obtaining-an-authentication-token-for-data-downloads).

    snakemake --configfile config.yml --use-conda

Before your first run, map your case IDs to file IDs using the `generate_sample_list` rule:

    snakemake --configfile config.yml --use-conda generate_sample_list

If you'd like only certain outputs, specify a target rule like so:

    snakemake --configfile config.yml --use-conda all_preprocess

Available target rules include: `all_download`, `all_preprocess`, `all_align`, `all_summary`

## About

fang-tcga is a backronym which stands for **f**inding **a**dditional **n**on-human **g**enomes in TCGA. It also uses code directly, as a submodule, or modified from, [hisss](https://github.com/louiejtaylor/hisss), a collab with ArwaAbbas.
