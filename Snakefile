# -*- mode: Snakemake -*-
#
# Workflow for dumping non-human reads from The Cancer Genome Atlas, mapping
# them to targets, then summarizing/visualizing results.
#

from pathlib import Path
import json, requests

def read_samples(sample_fp):
    with open(sample_fp, "r") as f:
        samples = f.read()
    return [sample.strip() for sample in samples.split("\n") if len(sample) > 0]

output_dir = Path(config["all"]["root_dir"])/Path(config["all"]["output_dir"])

try:
    Samples = read_samples(output_dir/"samples"/"sample_list.txt")
    print("Found",len(Samples),"samples")
except FileNotFoundError:
    print("Sample list not initialized: be sure to run with the generate_sample_list target rule before your first run.")

# Generate sample list (if necessary)
rule generate_sample_list:
    input:
        sample_file = config["download"]["acc_list_fp"]
    output:
        sample_list = output_dir/"samples"/"sample_list.txt"
    params:
        list_format=config["download"]["acc_list_format"],
        needs_mapping=config["download"]["map_sample_list_to_file_ids"],
        experimental_strategies=config["download"]["allowed_strategies"]
    run:
        if params.list_format=="list":
            acc_list = read_samples(input.sample_file)
        elif params.list_format=="gdc-json":
            json_list = json.load(open(input.sample_file))
            acc_list = [z["cases"][0]["case_id"] for z in json_list]
        else:
            raise Exception("Unknown sample list format: "+str(params.list_format))
        acc_list = list(set(acc_list))
        if params.needs_mapping:
            file_list = []
            for c in acc_list:
                json_mapping = json.loads(requests.get("https://api.gdc.cancer.gov/cases/"+c+"?expand=files").text)
                for fi in json_mapping["data"]["files"]:
                    if fi["data_format"]=="BAM":
                        if fi["experimental_strategy"] in params.experimental_strategies:
                            file_list.append([c, fi["file_id"], fi["file_name"]])
            with open(output_dir/"samples"/"sample_mapping.tsv", 'w') as o:
                o.write('\n'.join(['\t'.join(f) for f in file_list]))
            final_list = [f[1] for f in file_list]
        else:
            final_list = acc_list
        with open(output.sample_list, "w") as o:
            o.write('\n'.join(final_list))

# Acquire data
rule download_bams:
    input:
        sample_list = output_dir/"samples"/"sample_list.txt"
    output:
        output_dir/"download"/"{sample}.bam"
    params:
        token_str = ("-t "+ config["download"]["token_file"])*(len(config["download"]["token_file"]) > 0),
        sample = "{sample}",
        wdir = str(output_dir/"download"/"{sample}")
    conda:
        "download_env.yml"
    shell:
        """
        mkdir -p {params.wdir}
        gdc-client download -d {params.wdir} {params.token_str} {params.sample}
        mv {params.wdir}/https\:/api.gdc.cancer.gov/data/{params.sample}/*.bam {output}
        """

rule all_download:
    input:
        expand(str(output_dir/"download"/"{sample}.bam"), sample = Samples)
        
