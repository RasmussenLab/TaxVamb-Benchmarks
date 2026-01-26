import pandas as pd
import collections
import os
from pathlib import Path
import sys

# Set the directory the snakefile exists in. This makes us able to call the pipeline with the relevant src files from other directories.
THIS_FILE_DIR = Path(workflow.basedir)
# Define the src directory for the files used in the snakemake workflow
SRC_DIR = THIS_FILE_DIR / "files_used_in_snakemake_workflow"

# If the configfile is not set explicit fall back to the default
CONFIG_PATH = THIS_FILE_DIR / "config/config.yaml"
# TODO should go to configfile
if CONFIG_PATH.exists():
    configfile: CONFIG_PATH

# Get the output_directory defined by the user or fallback to current directory, which is the default way snakemake handles output directories
OUTDIR = Path("") if config.get("output_directory") is None else Path(config.get("output_directory"))

#### Setting parameters from the config file ####
##  For a more throughout description of what the different config options mean see the /config/config.yaml file

# Default resources used
default_walltime = config.get("default_walltime", "48:00:00")
default_threads = config.get("default_threads", 16)
print(default_threads)
default_mem_gb = config.get("default_mem_gb", 50)
default_gpu = ""

# Functions to get the config-defined threads/walltime/mem_gb for a rule and if not defined the default
threads_fn = lambda rulename: config.get(rulename, {"threads": default_threads}).get("threads", default_threads)
walltime_fn  = lambda rulename: config.get(rulename, {"walltime": default_walltime}).get("walltime", default_walltime)
mem_gb_fn  = lambda rulename: config.get(rulename, {"mem_gb": default_mem_gb}).get("mem_gb", default_mem_gb)
gpu_fn  = lambda rulename: config.get(rulename, {"gpu": default_gpu}).get("gpu", default_gpu)

# Assert that input files are actually passed to snakemake
if config.get("read_assembly_dir") == None:
    print("ERROR: read_assembly_dir not passed to snakemake as config. Define either. Eg. snakemake <arguments> --config read_file=<read file>. If in doubt refer to the README.md file")
    sys.exit()

# Set default values for dictonaries containg information about the input information
# The way snakemake parses snakefiles means we have to define them even though they will always be present
sample_id = dict()
sample_id_path= dict()
sample_id_path_assembly = dict()
sample_id_contig = collections.defaultdict()

read_fw = ""
read_rv = ""

df = pd.read_csv(config["read_assembly_dir"], sep=r"\s+", comment="#")
sample_id = collections.defaultdict(list)
sample_id_path = collections.defaultdict(dict)
sample_id_path_assembly = collections.defaultdict(dict)
for id, (sample, read1, read2, contig) in enumerate(zip( df["sample"], df.read1, df.read2, df.contig)):
    id = f"{Path(read1).stem.split("_")[0]}" # TODO: Change to be more robust
    sample_id[sample].append(id)
    sample_id_path[sample][id] = [read1, read2]
    sample_id_path_assembly[sample][id] = [contig]


# Setting the output paths for the user defined SPades files
contigs =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0])
read_fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
read_rv =  lambda wildcards: sample_id_path[wildcards.key][wildcards.id][1]

bamfiles_before =  str(OUTDIR) +  "/{key}/assembly_mapping_output/mapped/{id}.bam"

INDEX_SIZE = "12G"
MIN_CONTIG_LEN = 2000

print(sample_id)
rule all:
    input:
        [OUTDIR / bamfiles_before.format(key=k, id=i) for k, v in sample_id.items() for i in v]
        # expand(bamfiles_before, key=sample_id.keys(), id=sample_id.values())
        # expand(OUTDIR / "{key}/assembly_mapping_output/contigs.flt.fna.gz", key=sample_id.keys())


# Rename the contigs to keep sample information for later use 
rulename = "rename_contigs"
rule rename_contigs:
    input:
        contigs,
    output:
        OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.renamed.fasta"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    shell:
        """
        if [[ "{input}" == *.gz ]]; then
            zcat {input} | sed 's/^>/>S{wildcards.id}C/' > {output} 2> {log}
        else
            sed 's/^>/>S{wildcards.id}C/' {input} > {output} 2> {log}
        fi
        """

rulename = "cat_contigs"
rule cat_contigs:
    input: 
        lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.renamed.fasta", key=wildcards.key, id=sample_id[wildcards.key]),

    output: OUTDIR / "{key}/assembly_mapping_output/contigs.flt.fna.gz"
    threads: threads_fn(rulename)
    params: script =  SRC_DIR / "concatenate.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell: 
        "python {params.script} {output} {input} --keepnames -m {MIN_CONTIG_LEN} &> {log} "

# Index resulting contig-file with minimap2
rulename = "index"
rule index:
    input:
        contigs = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
    output:
        mmi = os.path.join(OUTDIR, "{key}", "contigs.flt.mmi")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: 
        THIS_FILE_DIR / "envs/minimap2.yaml"
    shell:
        "minimap2 -I {INDEX_SIZE} -d {output} {input} 2> {log}"

# This rule creates a SAM header from a FASTA file.
# We need it because minimap2 for truly unknowable reasons will write
# SAM headers INTERSPERSED in the output SAM file, making it unparseable.
# To work around this mind-boggling bug, we remove all header lines from
# minimap2's SAM output by grepping, then re-add the header created in this
# rule.
rulename = "dict"
rule dict:
    input:
        contigs = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
    output:
        dict = os.path.join(OUTDIR,"{key}", "contigs.flt.dict")  
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}" + rulename
    conda:
        THIS_FILE_DIR / "envs/samtools.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

# Generate bam files 
rulename = "minimap"
rule minimap:
    input:
        fw = read_fw,
        rv = read_rv,
        mmi = os.path.join(OUTDIR, "{key}", "contigs.flt.mmi"),
        dict = os.path.join(OUTDIR,"{key}", "contigs.flt.dict"),
    output:
        bam = bamfiles_before, 
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    conda:
        THIS_FILE_DIR / "envs/minimap2.yaml"
    shell:
        # See comment over rule "dict" to understand what happens here
        "minimap2 -t {threads} -ax sr {input.mmi} {input.fw} {input.rv} -N 5"
        " | grep -v '^@'"
        " | cat {input.dict} - "
        " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
        " > {output.bam} 2> {log}"


