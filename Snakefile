# -*- mode: Snakemake -*-
#
# Workflow for dumping non-human reads from The Cancer Genome Atlas, mapping
# them to targets, then summarizing/visualizing results.
#

from pathlib import Path

def read_samples(sample_fp):
    with open(sample_fp) as f:
        samples = f.read()
    return [sample.strip() for sample in samples.split(\n)]

Samples = read_samples(config["all"]["sample_key_fp"])
output_dir = Path(config["all"]["root_dir"])/Path(config["all"]["output_dir"])

# Acquire data
rule download_bams:
    output:
        output_dir/"download"/"{sample}.bam"
    shell:
        """
        # wget? curl? aspera? other NCI-specific tool?
        """

rule all_download:
    input:
        expand(str(output_dir/"download"/"{sample}.bam"), sample = Samples)
        
# Preprocess data
rule dump_unmapped_reads:
    input:
        rules.download_bams.output
    output:
        output_dir/"preproc"/"{sample}.fastq.gz"
    conda:
        "samtools_env.yml"
    shell:
        """
        
        """

rule all_dump_unmapped:
    input:
        expand(str(output_dir/"preproc"/"{sample}.fastq.gz"), sample = Samples)

# Align reads
rule build_aln_db:
    input:
        config["align"]["target_fasta"]
    output:
        # db for the aligner I use
    conda:
        "align_env.yml"
    shell:
        """

        """

rule align_reads:
    input:
        reads = rules.dump_unmapped_reads.output,
        db = rules.build_aln_db.output
    output:
        output_dir/"alignments"/"{sample}.bam"
    conda:
        "align_env.yml"
    shell:
        """

        """

rule all_align:
    input:
        expand(str(output_dir/"alignments"/"{sample}.bam"), sample = Samples)

# Summary
rule summarize_alignments:
    input:
        rules.align_reads.output
    output:
        # some kind of human and machine-readable summary--text file?
    conda:
        "samtools_env.yml"
    shell:
        """

        """

rule visualize_summary_data:
    input:
        # expand() statement for whatever `summarize_alignments` produces
    output:
        output_dir/"summary"/"all_summarized.Rmd"
    conda:
        "viz_env.yml"
    shell:
        """
        
        """

rule all_summarize:
    input:
        a = rules.visualize_summary_data.output,
        b = expand() #whatever rules.summarize_alignments produces
