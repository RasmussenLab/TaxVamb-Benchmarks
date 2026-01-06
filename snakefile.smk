import pandas as pd
import collections
import os
from pathlib import Path
import sys

# Set the directory the snakefile exists in. This makes us able to call the pipeline with the relevant src files from other directories.
THIS_FILE_DIR = Path(workflow.basedir)

# If the configfile is not set explicit fall back to the default
CONFIG_PATH = THIS_FILE_DIR / "config/config.yaml"
# TODO should go to configfile
if CONFIG_PATH.exists():
    configfile: CONFIG_PATH

# Define the src directory for the files used in the snakemake workflow
SRC_DIR = THIS_FILE_DIR / "files_used_in_snakemake_workflow"  

# Get the output_directory defined by the user or fallback to current directory, which is the default way snakemake handles output directories
OUTDIR = Path("") if config.get("output_directory") is None else Path(config.get("output_directory"))

#### Setting parameters from the config file ####
##  For a more throughout description of what the different config options mean see the /config/config.yaml file

# Default resources used
default_walltime = config.get("default_walltime", "48:00:00")
default_threads = config.get("default_threads", 16)
default_mem_gb = config.get("default_mem_gb", 50)
default_gpu = ""
use_minimap = config.get("use_minimap", True)
vamb_use_gpu = config.get("vamb_use_gpu")
vamb_extra_arg = ""
if vamb_use_gpu:
    vamb_extra_arg = "--cuda"
 
## ----------- ##

# Assert that input files are actually passed to snakemake
if config.get("bam_contig") == None:
    print("ERROR: bam_contig not passed to snakemake as config. Define either. Eg. snakemake <arguments> --config bam_contig=<bam contig>. If in doubt refer to the README.md file")
    sys.exit()

# Default paths for bamfiles 
bamfiles_before = OUTDIR / "{key}/assembly_mapping_output/mapped/{id}.bam"
bamfiles = OUTDIR / "{key}/assembly_mapping_output/mapped/{id}.bam"
contigs_all = OUTDIR / "{key}/assembly_mapping_output/contigs.flt.fna.gz"

# Set default values for dictonaries containg information about the input information
# The way snakemake parses snakefiles means we have to define them even though they will always be present
sample_id = dict()               
sample_id_path= dict() 
sample_id_path_assembly = dict()
sample_id_contig = collections.defaultdict()

# Read in the input data
if config.get("bam_contig") != None:
    df = pd.read_csv(config["bam_contig"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    contigs = collections.defaultdict(list)
    for id, (sample, bamfile, contig) in enumerate(zip(df["sample"], df.bamfile, df.contig)):
        id = f"sample{str(id)}"
        sample = sample
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [bamfile]
        sample_id_contig[sample] = contig
        contigs[sample].append(contig)
        bamfiles = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
        contigs_all = lambda wildcards: sample_id_contig[wildcards.key]

    # Check that all contigs are the same for each sample
    for sample in contigs.keys():
        for contig in contigs[sample]:
            if contigs[sample][0] != contig:
                print("Not all contigs are the same")
                sys.exit()

# Functions to get the config-defined threads/walltime/mem_gb for a rule and if not defined the default
threads_fn = lambda rulename: config.get(rulename, {"threads": default_threads}).get("threads", default_threads) 
walltime_fn  = lambda rulename: config.get(rulename, {"walltime": default_walltime}).get("walltime", default_walltime) 
mem_gb_fn  = lambda rulename: config.get(rulename, {"mem_gb": default_mem_gb}).get("mem_gb", default_mem_gb) 
gpu_fn  = lambda rulename: config.get(rulename, {"gpu": default_gpu}).get("gpu", default_gpu) 

# Sort the bam files 
rulename="sort"
rule sort:
    input:
        bamfiles,
    output:
        temp(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam"),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    conda: THIS_FILE_DIR / "envs/samtools.yaml"
    shell:
        """
	samtools sort --threads {threads} {input} -o {output} 2> {log}
	"""




rule all:
    input:
        # checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),







## Include the specific rules for each tool
include: THIS_FILE_DIR / "snakemake_modules/vamb_default.smk"
include: THIS_FILE_DIR / "snakemake_modules/metabat.smk"
include: THIS_FILE_DIR / "snakemake_modules/comebin.smk"
include: THIS_FILE_DIR / "snakemake_modules/metadecoder.smk"
include: THIS_FILE_DIR / "snakemake_modules/semibin.smk"
include: THIS_FILE_DIR / "snakemake_modules/taxvamb_using_mmseqs_classifications.smk"
include: THIS_FILE_DIR / "snakemake_modules/run_checkm_on_all_taxvamb.smk"
include: THIS_FILE_DIR / "snakemake_modules/gunc.smk"
# include: THIS_FILE_DIR / "snakemake_modules/va_extra_taxonomy_classifiers.smk"


## Collecting
all_bin_dirs_recluster = {
    # "kraken_taxvamb_default": OUTDIR / "{key}/kraken_taxvamb_default",
    # "run_taxvamb_centrifuge": OUTDIR / "{key}/centrifuge_taxvamb",
    # "run_taxvamb_gtdb": OUTDIR / "{key}/gtdb_taxvamb_default",
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown",
    # "metabuli_taxvamb_default": OUTDIR / "{key}/metabuli_taxvamb_default",
    # "kraken_taxvamb_no_predictor": OUTDIR / "{key}/kraken_taxvamb_no_predictor",
    # "centrifuge_taxvamb_no_predictor": OUTDIR / "{key}/centrifuge_taxvamb_no_predictor",
    # "run_taxvamb_gtdb_no_predictor": OUTDIR / "{key}/gtdb_taxvamb_default_no_predictor",
    # "metabuli_taxvamb_no_predictor": OUTDIR / "{key}/metabuli_taxvamb_no_predictor",
    # "default_vamb": OUTDIR / "{key}/vamb_default",
    # "kalmari_taxvamb_default": OUTDIR / "{key}/kalmari_taxvamb_default",
    # "trembl_taxvamb_default": OUTDIR / "{key}/trembl_taxvamb_default",
}
bin_dir_names_recluster = all_bin_dirs_recluster.keys()




        # reclustering = expand(OUTDIR / "{key}/reclustering/{bins_recluster}/output",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), 


        # checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), 

rule Checkm_reclustering:
    input:
        # With checkm2
        # directory = expand(OUTDIR / "{key}/checkm2/reclustering/{bins_recluster}",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), 
        checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),
# Run taxvamb 
rulename = "format_bins_class_recluster"
rule rename_vamb:
    input:
        bins = os.path.join(OUTDIR,"{key}",'vamb_default','vae_clusters_split.tsv'),
        latent = os.path.join(OUTDIR, '{key}','vamb_default/latent.npz'),
    output:
        bins = os.path.join(OUTDIR,"{key}",'vamb_default','vaevae_clusters_split.tsv'),
        latent = os.path.join(OUTDIR, '{key}','vamb_default/vaevae_latent.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    shell:
        """
            cp {input.bins} {output.bins}
            cp {input.latent} {output.latent}
        """

rulename = "recluster"
rule recluster:
    input: 
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = lambda wildcards: all_bin_dirs_recluster[wildcards.bins_recluster],
        bins = lambda wildcards: all_bin_dirs_recluster[wildcards.bins_recluster] / "vaevae_clusters_split.tsv",
        compo = lambda wildcards: all_bin_dirs_recluster[wildcards.bins_recluster] / "composition.npz",
    output:
        directory = directory(OUTDIR / "{key}/reclustering/{bins_recluster}/output"),
        headers = OUTDIR / "{key}/reclustering/{bins_recluster}/headers.txt",
        bins = OUTDIR / "{key}/reclustering/{bins_recluster}/output/clusters_reclustered.tsv",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    params:
        env = THIS_FILE_DIR / "reclustering"
    shell:
        """
        grep -E "^>" {input.contigs_decompressed} | cut -c 2- > {output.headers}
        rm -rf {output.directory}
        mkdir -p {output.directory}
        pixi run --manifest-path {params.env} start --num_process {threads} \
        {input.bins} {input.directory}/vaevae_latent.npz \
        {input.contigs_decompressed} \
        {output.headers} \
        {output.directory} \
        kmeans \
        """

# Run taxvamb 
rulename = "format_bins_class_recluster"
rule format_bins_class_recluster:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins = OUTDIR / "{key}/reclustering/{bins_recluster}/output/clusters_reclustered.tsv",
    output:
        directory = directory(OUTDIR / "{key}/reclustering/formatted_vamb_bins/{bins_recluster}"),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """


rulename = "checkm_class_recluster"
rule checkm_class_recluster:
    input: 
        bin_dir = directory(OUTDIR / "{key}/reclustering/formatted_vamb_bins/{bins_recluster}"),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/reclustering/{bins_recluster}"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """
