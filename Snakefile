# -*- mode: Snakemake -*-
#
# Workflow for dumping non-human reads from The Cancer Genome Atlas, mapping
# them to targets, then summarizing/visualizing results.
#

from pathlib import Path
import json, requests, os, time

def read_samples(sample_fp):
    with open(sample_fp, "r") as f:
        samples = f.read()
    return [sample.strip() for sample in samples.split("\n") if len(sample.strip()) > 0]

output_dir = Path(config["all"]["root_dir"])/Path(config["all"]["output_dir"])

try:
    Samples = read_samples(output_dir/"samples"/"sample_list.txt")
    print("Found",len(Samples),"samples")
except FileNotFoundError:
    print("Sample list not initialized: be sure to run with the generate_sample_list target rule before your first run.")
    Samples = []

print("Found " +str(len(Samples))+" samples.")

# Generate sample list (first run)
rule generate_sample_list:
    input:
        sample_file = config["download"]["acc_list_fp"]
    output:
        sample_list = output_dir/"samples"/"sample_list.txt"
    params:
        output_dir = output_dir/"samples"/"indiv",
        list_format=config["download"]["acc_list_format"],
        needs_mapping=config["download"]["map_sample_list_to_file_ids"],
        experimental_strategies=config["download"]["allowed_strategies"]
    run:
        shell("mkdir -p {params.output_dir}")

        # check format of sample list
        if params.list_format=="list":
            acc_list = read_samples(input.sample_file)
        elif params.list_format=="gdc-json":
            json_list = json.load(open(input.sample_file))
            acc_list = [z["cases"][0]["case_id"] for z in json_list]
        else: # handle gdc-manifest?
            raise Exception("Unknown sample list format: "+str(params.list_format))
        acc_list = list(set(acc_list))

        # only run if user has case ids and wants pipeline to find the associated files
        if params.needs_mapping:
            for c in acc_list:
                if os.path.isfile(params.output_dir/(str(c)+".txt")): # sample has been mapped
                    print("found map for "+c+", skipping")
                    continue

                # multiple tries often required for successful download. abstract to config?
                retries = 3
                while retries > 0:
                    try:
                        response = requests.get("https://api.gdc.cancer.gov/cases/"+c+"?expand=files", timeout = 10)
                        break
                    except requests.exceptions.Timeout:
                        print("connection failed for "+c+", retrying "+str(retries)+"x")
                        continue
                    except requests.exceptions.ConnectionError:
                        print("connectionError for "+c+", retrying "+str(retries)+"x")
                        time.sleep(3)
                        continue

                # find associated files in user-specified acceptable strategies
                json_mapping = json.loads(response.text)
                for fi in json_mapping["data"]["files"]:
                    if fi["data_format"]=="BAM":
                        if fi["experimental_strategy"] in params.experimental_strategies:
                            with open(params.output_dir/(str(c)+".txt"),"w") as o:
                                o.write(fi["file_id"]+"\n")
            # output sample list
            shell("cat {params.output_dir}/*.txt | sort | uniq > {output.sample_list}")
        else: # user supplied file ids, not case ids
            final_list = acc_list
            with open(output.sample_list, "w") as o:
                o.write('\n'.join(list(set(final_list))))

# Acquire data
rule download_bams:
    input:
        sample_list = output_dir/"samples"/"sample_list.txt"
    output:
        temp(output_dir/"download"/"{sample}.bam")
    params:
        token_str = ("-t "+ config["download"]["token_file"])*(len(config["download"]["token_file"]) > 0),
        sample = "{sample}",
        wdir = str(output_dir/"download"/"{sample}")
    conda:
        "envs/download_env.yml"
    threads: 1
    shell:
        """
        mkdir -p {params.wdir}
        gdc-client download -n {threads} -d {params.wdir} {params.token_str} {params.sample}
        mv {params.wdir}/https\:/api.gdc.cancer.gov/data/{params.sample}/*.bam {output}
        """

rule all_download:
    input:
        expand(str(output_dir/"download"/"{sample}.bam"), sample = Samples)
       
# Preprocess data
rule bam_to_reads:
    input:
        output_dir/"download"/"{sample}.bam"
    output:
        temp(output_dir/"reads"/"{sample}.fastq.gz")
    params:
        out_fastq = str(output_dir/"reads"/"{sample}.fastq")
    conda:
        "envs/alignment_processing_env.yml"
    shell:
        """
        bedtools bamtofastq -i {input} -fq {params.out_fastq}
        gzip {params.out_fastq}
        """

rule all_preprocess:
    input:
        expand(str(output_dir/"reads"/"{sample}.fastq.gz"), sample = Samples)

# Read alignments/db building
rule build_aln_db:
    input:
        config["align"]["target_fasta"]
    output:
        output_dir/"db"/"db.mmi"
    conda:
        "envs/align_env.yml"
    threads: 1
    shell:
        """
        minimap2 -x sr -t {threads} -d {output} {input}
        """

rule align_reads:
    input:
        reads = output_dir/"reads"/"{sample}.fastq.gz",
        db = output_dir/"db"/"db.mmi"
    output:
        temp(output_dir/"alignments"/"{sample}.sam")
    conda:
        "envs/align_env.yml"
    threads: 1
    shell:
        """
        minimap2 -t {threads} -ax sr {input.db} {input.reads} > {output}
        """

rule all_align:
    input:
        expand(str(output_dir/"alignments"/"{sample}.bam"), sample = Samples)

rule process_alignment:
    input:
        output_dir/"alignments"/"{sample}.sam"
    output:
        bam = temp(output_dir/"alignments"/"{sample}.bam"),
        sorted = output_dir/"alignments"/"{sample}.sorted.bam",
        bai = output_dir/"alignments"/"{sample}.sorted.bam.bai"
    conda:
        "envs/alignment_processing_env.yml"
    params:
        target = config["align"]["target_fasta"]
    threads: 1
    shell:
        """
        samtools view -@ {threads} -bT {params.target} {input} > {output.bam}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        samtools index -@ {threads} {output.sorted} {output.bai}
        """

# Summarize alignments
rule mapping_stats:
    input: 
        bam = output_dir/"alignments"/"{sample}.sorted.bam",
        idx = output_dir/"alignments"/"{sample}.sorted.bam.bai"
    output:
        output_dir/"summary"/"idxstats"/"{sample}.sorted.idxstats.tsv"
    conda:
        "envs/alignment_processing_env.yml"
    threads: 1
    shell:
        """
        samtools idxstats -@ {threads} {input.bam} > {output}
        """

rule calculate_coverage:
    input:
        bam = output_dir/"alignments"/"{sample}.sorted.bam",
        idx = output_dir/"alignments"/"{sample}.sorted.bam.bai"
    output:
        output_dir/"summary"/"coverage"/"{sample}.genomecoverage.txt"
    conda:
        "envs/alignment_processing_env.yml"
    threads: 1
    shell:
        """
        samtools view -@ {threads} -b {input.bam}|genomeCoverageBed -ibam stdin|grep -v 'genome'| perl hisss/scripts/coverage_counter.pl > {output}	
        """

rule combine_coverage_stats:
    input:
        cov = output_dir/"summary"/"coverage"/"{sample}.genomecoverage.txt",
        stats = output_dir/"summary"/"idxstats"/"{sample}.sorted.idxstats.tsv"
    output:
        output_dir/"summary"/"individual"/"{sample}.align.summary.txt"
    conda:
        "envs/r_plot.yml"
    shell:
        """
        Rscript hisss/scripts/summarize_alignments.R {input.cov} {input.stats} {output}
        """

#Combine information for all samples into a single file
rule all_summary:
    input:
        expand(str(output_dir/"summary"/"individual"/"{sample}.align.summary.txt"),sample=Samples)
    output:
        output_dir/"summary"/"all_align_summary.txt"
    params:
        align_dir = output_dir/"summary"/"individual"
    shell:
        """
        echo -e "Sample\tAlignTarget\tFractionCoverage\tTargetLength\tMappedReads" > {output}
        cat {params.align_dir}/*.align.summary.txt >> {output}
        """


rule plot_depth:
    input:
        output_dir/"alignments"/"{sample}.sorted.bam"
    output:
        cov=output_dir/"summary"/"plots"/"{sample}.cov.depth.txt"
    params:
        plot=output_dir/"summary"/"plots"/"{sample}.cov.pdf"
    conda:
        "envs/r_plot.yml"
    threads: 1
    shell:
        """
        samtools depth -@ {threads} -a {input} > {output.cov}	
        if [ -s {output.cov} ]; then
            Rscript hisss/scripts/plot_coverage.R {output.cov} {params.plot};
        else
            echo "No valid alignments detected";
        fi
        """

#Combine information for all samples into a single file
rule all_plot:
    input:
        expand(str(output_dir/"summary"/"plots"/"{sample}.cov.depth.txt"),sample=Samples)
    output:
        output_dir/"summary"/"all_plot_summary.txt"
    params:
        plot_dir = output_dir/"summary"/"plots"
    shell:
        """
        cat {params.plot_dir}/*.cov.depth.txt > {output}
        """

rule all:
    input:
        a=output_dir/"summary"/"all_align_summary.txt",
        b=output_dir/"summary"/"all_plot_summary.txt"
