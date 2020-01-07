# -*- mode: Snakemake -*-
#
# Workflow for dumping non-human reads from The Cancer Genome Atlas, mapping
# them to targets, then summarizing results.
#

def read_samples(sample_fp):
    with open(sample_fp) as f:
        samples = f.read()
    return [sample.strip() for sample in samples.split(\n)]

Samples = read_samples(config["sample_key_fp"])

rule download_bams:

rule dump_unmapped_reads:

rule bam_to_fastq:

rule all_dump_unmapped:

rule build_align_db:

rule align_reads:

rule all_align:

rule summarize_alignments:

rule visualize_summary_data:

rule all:
