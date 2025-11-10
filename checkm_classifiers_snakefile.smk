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

# Minimum contig length used
MIN_CONTIG_LEN = 99999999999999# int(config.get("min_contig_len", 2000)) 

# Other options
CUDA = True if config.get("cuda") ==  "True" else False
 
## ----------- ##

# Assert that input files are actually passed to snakemake
if config.get("read_file") == None and config.get("read_assembly_dir") == None and config.get("bam_contig") == None:
    print("ERROR: read_file or read_assembly_dir or bam_contig not passed to snakemake as config. Define either. Eg. snakemake <arguments> --config read_file=<read file>. If in doubt refer to the README.md file")
    sys.exit()

# Set default paths for the SPades outputfiles - running the pipeline from allready assembled reads overwrite these values
contigs =  OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.fasta"
contigs_paths =  OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.paths"
assembly_graph = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/assembly_graph_after_simplification.gfa"

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

read_fw = ""
read_rv = ""
# If the read_file is defined the pipeline will also run SPades and assemble the reads
if config.get("read_file") != None:
    df = pd.read_csv(config["read_file"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    for id, (sample, read1, read2) in enumerate(zip(df["sample"], df.read1, df.read2)):
        id = f"sample{str(id)}"
        sample = sample
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [read1, read2]
    read_fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
    read_rv =  lambda wildcards: sample_id_path[wildcards.key][wildcards.id][1]

# If read_assembly dir is defined the pipeline will run user defined SPades output files
if config.get("read_assembly_dir") != None:
    df = pd.read_csv(config["read_assembly_dir"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    sample_id_path_assembly = collections.defaultdict(dict)
    for id, (sample, read1, read2, assembly) in enumerate(zip( df["sample"], df.read1, df.read2, df.assembly_dir)):
        id = f"sample{str(id)}"
        sample = sample
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [read1, read2]
        sample_id_path_assembly[sample][id] = [assembly]

    # Setting the output paths for the user defined SPades files
    contigs =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0]) / "contigs.fasta"
    assembly_graph  =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0]) / "assembly_graph_after_simplification.gfa"
    contigs_paths  =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0]) / "contigs.paths"

    read_fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
    read_rv =  lambda wildcards: sample_id_path[wildcards.key][wildcards.id][1]

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

print(sample_id_contig)
            


# Functions to get the config-defined threads/walltime/mem_gb for a rule and if not defined the default
threads_fn = lambda rulename: config.get(rulename, {"threads": default_threads}).get("threads", default_threads) 
walltime_fn  = lambda rulename: config.get(rulename, {"walltime": default_walltime}).get("walltime", default_walltime) 
mem_gb_fn  = lambda rulename: config.get(rulename, {"mem_gb": default_mem_gb}).get("mem_gb", default_mem_gb) 
gpu_fn  = lambda rulename: config.get(rulename, {"gpu": default_gpu}).get("gpu", default_gpu) 

rulename = "all"
rule all:
    input:
        # Biodata (3)
        checkm = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),
        gunc = expand(OUTDIR /  "{key}/tmp/gunc.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        gunc_2 = expand(OUTDIR /  "{key}/tmp/gunc_2.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        checkm_kalmari = expand(OUTDIR /  "{key}/checkm2/kalmari_taxvamb",key=sample_id.keys()), # Only for the 3 new biological samples
        checkm_trembl = expand(OUTDIR /  "{key}/checkm2/trembl_taxvamb",key=sample_id.keys()),# Only for the 3 new biological samples
        # running
        checkm_gtdb = expand(OUTDIR /  "{key}/checkm2/gtdb_taxvamb",key=sample_id.keys()),
        # from here to run
        # checkm_centrifuge = expand(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb_no_predictor",key=sample_id.keys()), 
        # checkm_kraken = expand(OUTDIR /  "{key}/checkm2/kraken_taxvamb_no_predictor",key=sample_id.keys()),
        # checkm_metabuli = expand(OUTDIR /  "{key}/checkm2/metabuli_taxvamb_no_predictor",key=sample_id.keys()), 
        # checkm_centrifuge = expand(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb",key=sample_id.keys()), 
        # checkm_kraken = expand(OUTDIR /  "{key}/checkm2/kraken_taxvamb",key=sample_id.keys()),
        # checkm_metabuli = expand(OUTDIR /  "{key}/checkm2/metabuli_taxvamb",key=sample_id.keys()), 
        # metabuli= expand(OUTDIR /  "{key}/classifiers/metabuli",key=sample_id.keys()),
        # kraken2 = expand(OUTDIR /  "{key}/classifiers/kraken2/kraken2_predictions.tsv",key=sample_id.keys()),
        # centrifuge = expand(OUTDIR /  "{key}/classifiers/centrifuge/centrifuge_predictions.tsv",key=sample_id.keys()),
        # checkm_gtdb = expand(OUTDIR /  "{key}/checkm2/gtdb_taxvamb",key=sample_id.keys()),
        # bins_gtdb =     expand(os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default','vaevae_clusters_unsplit.tsv'),key=sample_id.keys()),
        # bins_kalmari =  expand(os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default','vaevae_clusters_unsplit.tsv'),key=sample_id.keys()),
        # bins_trembl =   expand(os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default','vaevae_clusters_unsplit.tsv'),key=sample_id.keys()),
        # bins_trembl_fasta =   expand(os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default_bins'),key=sample_id.keys()),
        # bins_gtdb_fasta =   expand(os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default_bins'),key=sample_id.keys()),
        # bins_kalmari_fasta =   expand(os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default_bins'),key=sample_id.keys()),
        # checkm_gtdb = expand(OUTDIR /  "{key}/checkm2/gtdb_taxvamb",key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        # comebin = expand(OUTDIR /  "{key}/comebin", key=sample_id.keys()),
        # vamb_default = expand( OUTDIR / "{key}" / 'vamb_default' / 'vae_clusters_unsplit.tsv', key=sample_id.keys()),
        # metadecoder = expand(OUTDIR / "{key}/metadecoder/seed.seed", key=sample_id.keys()),
        # metabat = expand(OUTDIR /  "{key}/metabat/depht.txt", key=sample_id.keys()),
        # semibin = expand(OUTDIR /  "{key}/semibin", key=sample_id.keys()),
        # metabuli=   expand( OUTDIR /  "{key}/classifiers/metabuli",    key=sample_id.keys()),
        # mmseqs2 = expand(OUTDIR /  "{key}/classifiers/mmseqs2", key=sample_id.keys()),
        # kraken2 =  expand( OUTDIR /  "{key}/classifiers/kraken2",          key=sample_id.keys()),     
        # centrifuge = expand(OUTDIR /  "{key}/classifiers/centrifuge",           key=sample_id.keys()),




output_name = "gt_tax/{gt_tax}"
all_targets = [
     # "no_unknown_old_version_taxvamb_formatted_classifications.tsv"
    # "gt_tax_airways_removed_10pct.tsv",
    # "gt_tax_airways_removed_25pct.tsv",
    # "gt_tax_airways_removed_50pct.tsv",
    # "gt_tax_airways.tsv",
    # "gt_tax_airways_unknown_10pct.tsv",
    # "gt_tax_airways_unknown_25pct.tsv",
    # "gt_tax_airways_unknown_50pct.tsv",
    #"less_than_q_0.0_centrifuge_black_sea_clas.tsv"   ,
    #"less_than_q_100.0_centrifuge_black_sea_clas.tsv",
    #"less_than_q_144.0_centrifuge_black_sea_clas.tsv",
    #"less_than_q_196.0_centrifuge_black_sea_clas.tsv",
    #"less_than_q_324.0_centrifuge_black_sea_clas.tsv",
    #"less_than_q_64.0_centrifuge_black_sea_clas.tsv",
    #"less_than_q_676.0_centrifuge_black_sea_clas.tsv",
    #"less_than_q_81.0_centrifuge_black_sea_clas.tsv",
    "gt_tax_airways_shuffled_0pct.tsv",
    "gt_tax_airways_shuffled_100pct.tsv",
    "gt_tax_airways_shuffled_10pct.tsv",
    "gt_tax_airways_shuffled_25pct.tsv",
    "gt_tax_airways_shuffled_50pct.tsv",
    "gt_tax_airways_shuffled_75pct.tsv",
]




rule Airways_test_different_classifications:
    input:
        test = expand(OUTDIR /  "{key}" / output_name / "vaevae_clusters_split.tsv",gt_tax = all_targets, key=sample_id.keys()),
        # outdir = expand(OUTDIR /  "{key}/checkm2" / output_name, gt_tax = all_targets, key=sample_id.keys()),

# Run taxvamb 
rulename = "run_taxvamb_kraken"
rule run_taxvamb_differenct_classifiers:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "gt_tax/{gt_tax}",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", output_name)),
        bins = os.path.join(OUTDIR,"{key}",output_name,'vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}',output_name,'composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} 
        """

# Run taxvamb 
rulename = "format_bins_class"
rule format_bins_class_2:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins =  os.path.join(OUTDIR,"{key}",output_name,'vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", "format_clusters",output_name)),
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

rulename = "checkm_class"
rule checkm_class_2:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", "format_clusters",output_name)),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2" / output_name),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """



rule Biotest_Airways:
    input:
        # checkm = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        test = expand(OUTDIR /  "{key}/metabuli_taxvamb_default/vaevae_clusters_split.tsv",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        #test2 = expand(OUTDIR /  "{key}/1_metabuli_taxvamb_default/vaevae_clusters_split.tsv",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING

rule Classify_genomad:
    input:
        # expand(os.path.join(OUTDIR,"{key}",'geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv'), key=sample_id.keys())
        outdir = expand(OUTDIR /  "{key}/checkm2/plasmid_vira_split", key=sample_id.keys())

rule Classify_genomad_2:
    input:
        expand(os.path.join(OUTDIR,"{key}",'geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv'), key=sample_id.keys())


rule mmseqs_special:
    input:
        out1 = expand(OUTDIR /  "{key}/checkm2/kalmari_taxvamb", key=sample_id.keys()),
        out2 = expand(OUTDIR /  "{key}/checkm2/trembl_taxvamb", key=sample_id.keys())

geNomad_db = THIS_FILE_DIR / "genomad_db" / "genomad_db"

rulename = "download_genomad_db"
rule download_genomad_db:
    output:
        geNomad_db_path = directory(geNomad_db),
        dir_geNomad_db = directory(THIS_FILE_DIR / "genomad_db"),
        geNomad_db = protected(THIS_FILE_DIR / "genomad_db" / "genomad_db" / "genomad_db"),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/genomad.yaml"
    shell:
        """
        mkdir -p {output.dir_geNomad_db}
        genomad download-database {output.dir_geNomad_db}
        """

# 9. Run geNomad to get contig plasmid, chromosome, and virus classificaiton scores
rulename = "run_geNomad"
rule run_geNomad:
    input:
        contigs = contigs_all,
        geNomad_db = geNomad_db
    output:
        directory(os.path.join(OUTDIR,"{key}",'geNomad')),
        os.path.join(OUTDIR,"{key}",'rule_completed_checks/run_geNomad.finished'),
        os.path.join(OUTDIR,"{key}",'geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv')
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/genomad.yaml"
    shell:
        """
        # genomad end-to-end --cleanup {input.contigs} {output[0]} {input.geNomad_db} --threads {threads} &> {log}
        genomad end-to-end {input.contigs} {output[0]} {input.geNomad_db} --threads {threads} &> {log}
        touch {output[1]}
        """


# Run taxvamb 
rulename = "format_bins_class"
rule format_bins_class_3:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/transfer/plasmid_vira_split.tsv"
    output:
        directory = directory(OUTDIR / "{key}/plasmid_vira_split_formatted_vamb_bins"),
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


rulename = "checkm_class"
rule checkm_class_3:
    input: 
        bin_dir = directory(OUTDIR / "{key}/plasmid_vira_split_formatted_vamb_bins"),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/plasmid_vira_split"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """


rule Biotest:
    input:
        # checkm = expand(OUTDIR /  "{key}/tmp/run_all.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING

## Collecting
all_bin_dirs_recluster = {
    "kraken_taxvamb_default": OUTDIR / "{key}/kraken_taxvamb_default",
    "run_taxvamb_centrifuge": OUTDIR / "{key}/centrifuge_taxvamb",
    # "run_taxvamb_gtdb": OUTDIR / "{key}/gtdb_taxvamb_default",
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown",
    "metabuli_taxvamb_default": OUTDIR / "{key}/metabuli_taxvamb_default",
    # "kraken_taxvamb_no_predictor": OUTDIR / "{key}/kraken_taxvamb_no_predictor",
    # "centrifuge_taxvamb_no_predictor": OUTDIR / "{key}/centrifuge_taxvamb_no_predictor",
    # "run_taxvamb_gtdb_no_predictor": OUTDIR / "{key}/gtdb_taxvamb_default_no_predictor",
    # "metabuli_taxvamb_no_predictor": OUTDIR / "{key}/metabuli_taxvamb_no_predictor",
    "default_vamb": OUTDIR / "{key}/vamb_default",
    "kalmari_taxvamb_default": OUTDIR / "{key}/kalmari_taxvamb_default",
    "trembl_taxvamb_default": OUTDIR / "{key}/trembl_taxvamb_default",
}
bin_dir_names_recluster = all_bin_dirs_recluster.keys()

rule mmseqs_with_gunc:
    input:
        # 
        # For CAMI without checkm
        directory = expand(OUTDIR / "{key}/reclustering/{bins_recluster}/output",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # With checkm
        # 
        ## !! These are used on jun 14 !!
        # checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # directory = expand(OUTDIR / "{key}/checkm2/reclustering/{bins_recluster}",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), 
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),
        # gunc = expand(OUTDIR /  "{key}/tmp/gunc.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING

rule run_comebin_checkm:
    input:
        checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),

rule metabat_cami:
    input:
        checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),

rule comebin_madnessReclustering:
    input:
        checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # bins = expand(OUTDIR / "{key}/comebin/comebin_res/comebin_res_bins", key=sample_id.keys())

rule metadecoder_madnessReclustering:
    input:
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),

rule madnessReclustering:
    input:
        checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        # 
        # For CAMI without checkm
        # directory = expand(OUTDIR / "{key}/reclustering/{bins_recluster}/output",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # With checkm
        # 
        ## !! These are used on jun 14 !!
        # checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # directory = expand(OUTDIR / "{key}/checkm2/reclustering/{bins_recluster}",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), 
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),

rule aeroReclustering:
    input:
        checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        # 
        # For CAMI without checkm
        # directory = expand(OUTDIR / "{key}/reclustering/{bins_recluster}/output",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # With checkm
        # 
        ## !! These are used on jun 14 !!
        # checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # directory = expand(OUTDIR / "{key}/checkm2/reclustering/{bins_recluster}",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), 
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),

rule Reclustering:
    input:
        # checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # For CAMI without checkm
        # directory = expand(OUTDIR / "{key}/reclustering/{bins_recluster}/output",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # With checkm
        # 
        ## !! These are used on jun 14 !!
        checkm2 = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # directory = expand(OUTDIR / "{key}/checkm2/reclustering/{bins_recluster}",key=sample_id.keys(), bins_recluster=bin_dir_names_recluster), 
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),


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
        # hmm_path = OUTDIR / "marker_data/{key}/markers.hmmout",
    output:
        directory = directory(OUTDIR / "{key}/reclustering/{bins_recluster}/output"),
        headers = OUTDIR / "{key}/reclustering/{bins_recluster}/headers.txt",
        bins = OUTDIR / "{key}/reclustering/{bins_recluster}/output/clusters_reclustered.tsv",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    params:
        env = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/pixi_vamb/misc_scripts/reclustering"
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
        # --hmmout_path \
        # {input.hmm_path} \

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

# rule test:
#     input: 
#         contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
#         directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
#         bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_default','vaevae_clusters_split.tsv'),
#         compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_default/composition.npz'),
#         hmm_path = OUTDIR / "marker_data/{key}/markers.hmmout"
#     output:
#         directory = directory(OUTDIR / "{key}/reclustering/metabuli/output"),
#         headers = OUTDIR / "{key}/reclustering/headers.txt"
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     params:
#         env = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/pixi_vamb/misc_scripts/reclustering"
#     shell:
#         """
#         grep -E "^>" {input.contigs_decompressed} | cut -c 2- > {output.headers}
#         mkdir -p {output.directory}
#         pixi run --manifest-path {params.env} start --num_process {threads} \
#         --hmmout_path \
#         {input.hmm_path} \
#         {input.bins} {input.directory}/vaevae_latent.npz \
#         {input.contigs_decompressed} \
#         {output.headers} \
#         {output.directory} \
#         kmeans \
#         """

# rule test:
#     input: 
#         contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
#         directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
#         bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_default','vaevae_clusters_split.tsv'),
#         compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_default/composition.npz'),
#         hmm_path = OUTDIR / "marker_data/{key}/markers.hmmout"
#     output:
#         directory = directory(OUTDIR / "{key}/reclustering/output"),
#         headers = OUTDIR / "{key}/reclustering/headers.txt"
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     params:
#         env = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/pixi_vamb/misc_scripts/reclustering"
#     shell:
#         """
#         grep -E "^>" {input.contigs_decompressed} | cut -c 2- > {output.headers}
#         mkdir -p {output.directory}
#         pixi run --manifest-path {params.env} start --num_process {threads} \
#         {input.bins} {input.directory}/vaevae_latent.npz \
#         {input.contigs_decompressed} \
#         {output.headers} \
#         {output.directory} \
#         kmeans \
#         """
#         # --hmmout_path \
#         # {input.hmm_path} \

rule airways:
    input:
        # Biodata (3)
        # checkm = expand(OUTDIR /  "{key}/tmp/run_all.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        test = expand(OUTDIR /  "{key}/metabuli_taxvamb_default/vaevae_clusters_split.tsv",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        test2 = expand(OUTDIR /  "{key}/kraken_taxvamb_default/vaevae_clusters_split.tsv",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # test3 = expand(OUTDIR /  "{key}/gtdb_taxvamb_default/vaevae_clusters_split.tsv",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        test4 = expand(OUTDIR /  "{key}/centrifuge_taxvamb/vaevae_clusters_split.tsv",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),
        # gunc = expand(OUTDIR /  "{key}/tmp/gunc.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # gunc_2 = expand(OUTDIR /  "{key}/tmp/gunc_2.done",key=sample_id.keys()), ## WHEN USING THIS MAKE SURE NOT TO RERUN ANYTHING
        # checkm_kalmari = expand(OUTDIR /  "{key}/checkm2/kalmari_taxvamb",key=sample_id.keys()), # Only for the 3 new biological samples
        # checkm_trembl = expand(OUTDIR /  "{key}/checkm2/trembl_taxvamb",key=sample_id.keys()),# Only for the 3 new biological samples
        # running
        # checkm_gtdb = expand(OUTDIR /  "{key}/checkm2/gtdb_taxvamb",key=sample_id.keys()),
        # from here to run
        # checkm_centrifuge = expand(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb_no_predictor",key=sample_id.keys()), 
        # checkm_kraken = expand(OUTDIR /  "{key}/checkm2/kraken_taxvamb_no_predictor",key=sample_id.keys()),
        # checkm_metabuli = expand(OUTDIR /  "{key}/checkm2/metabuli_taxvamb_no_predictor",key=sample_id.keys()), 
        # checkm_centrifuge = expand(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb",key=sample_id.keys()), 
        # checkm_kraken = expand(OUTDIR /  "{key}/checkm2/kraken_taxvamb",key=sample_id.keys()),
        # checkm_metabuli = expand(OUTDIR /  "{key}/checkm2/metabuli_taxvamb",key=sample_id.keys()), 
        # metabuli= expand(OUTDIR /  "{key}/classifiers/metabuli",key=sample_id.keys()),
        # kraken2 = expand(OUTDIR /  "{key}/classifiers/kraken2/kraken2_predictions.tsv",key=sample_id.keys()),
        # centrifuge = expand(OUTDIR /  "{key}/classifiers/centrifuge/centrifuge_predictions.tsv",key=sample_id.keys()),
        # checkm_gtdb = expand(OUTDIR /  "{key}/checkm2/gtdb_taxvamb",key=sample_id.keys()),
        # bins_gtdb =     expand(os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default','vaevae_clusters_unsplit.tsv'),key=sample_id.keys()),
        # bins_kalmari =  expand(os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default','vaevae_clusters_unsplit.tsv'),key=sample_id.keys()),
        # bins_trembl =   expand(os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default','vaevae_clusters_unsplit.tsv'),key=sample_id.keys()),
        # bins_trembl_fasta =   expand(os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default_bins'),key=sample_id.keys()),
        # bins_gtdb_fasta =   expand(os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default_bins'),key=sample_id.keys()),
        # bins_kalmari_fasta =   expand(os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default_bins'),key=sample_id.keys()),
        # checkm_gtdb = expand(OUTDIR /  "{key}/checkm2/gtdb_taxvamb",key=sample_id.keys()),
        # checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        # checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        # checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        # checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        # comebin = expand(OUTDIR /  "{key}/comebin", key=sample_id.keys()),
        # vamb_default = expand( OUTDIR / "{key}" / 'vamb_default' / 'vae_clusters_unsplit.tsv', key=sample_id.keys()),
        # metadecoder = expand(OUTDIR / "{key}/metadecoder/seed.seed", key=sample_id.keys()),
        # metabat = expand(OUTDIR /  "{key}/metabat/depht.txt", key=sample_id.keys()),
        # semibin = expand(OUTDIR /  "{key}/semibin", key=sample_id.keys()),
        # metabuli=   expand( OUTDIR /  "{key}/classifiers/metabuli",    key=sample_id.keys()),
        # mmseqs2 = expand(OUTDIR /  "{key}/classifiers/mmseqs2", key=sample_id.keys()),
        # kraken2 =  expand( OUTDIR /  "{key}/classifiers/kraken2",          key=sample_id.keys()),     
        # centrifuge = expand(OUTDIR /  "{key}/classifiers/centrifuge",           key=sample_id.keys()),

#### Rules general for all tools ####

# If only reads are passed run metaspades to assemble the reads
rulename = "spades"
rule spades:
    input:
       fw = read_fw, 
       rv = read_rv, 
    output:
       outdir = directory(OUTDIR / "{key}/assembly_mapping_output/spades_{id}"),
       outfile = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.fasta",
       graph = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/assembly_graph_after_simplification.gfa",
       graphinfo  = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.paths",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    conda: THIS_FILE_DIR / "envs/spades_env.yaml"
    shell:
       "spades.py --meta "
       "-t {threads} -m 180 "
       "-o {output.outdir} -1 {input.fw} -2 {input.rv} " 
       "-t {threads} --memory {resources.mem_gb} &> {log} " 

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
        sed 's/^>/>S{wildcards.id}C/' {input} > {output} 2> {log}
        """

# Cat the contigs together in one file to later map each pair of reads against all the contigs together
rulename="cat_contigs"
rule cat_contigs:
    input: lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.renamed.fasta", key=wildcards.key, id=sample_id[wildcards.key]),
    output: OUTDIR / "{key}/assembly_mapping_output/contigs.flt.fna.gz"
    threads: threads_fn(rulename)
    params: script =  SRC_DIR / "concatenate.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell: 
        "python {params.script} {output} {input} --keepnames -m {MIN_CONTIG_LEN} &> {log} "  


# if use_minimap:
#     INDEX_SIZE = "12G" # get_config("index_size", "12G", r"[1-9]\d*[GM]$")
#     # Index resulting contig-file with minimap2
#     rule index:
#         input:
#             # contigs = os.path.join(OUTDIR,"contigs.flt.fna.gz")
#             contigs = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
#         output:
#             mmi = os.path.join(OUTDIR, "{key}", "contigs.flt.mmi")
#         threads: threads_fn(rulename)
#         resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#         benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#         log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#         conda: 
#             THIS_FILE_DIR / "envs/minimap2.yaml"
#         shell:
#             "minimap2 -I {INDEX_SIZE} -d {output} {input} 2> {log}"
#
#     # This rule creates a SAM header from a FASTA file.
#     # We need it because minimap2 for truly unknowable reasons will write
#     # SAM headers INTERSPERSED in the output SAM file, making it unparseable.
#     # To work around this mind-boggling bug, we remove all header lines from
#     # minimap2's SAM output by grepping, then re-add the header created in this
#     # rule.
#     rule dict:
#         input:
#             # contigs = os.path.join(OUTDIR,"contigs.flt.fna.gz")
#             contigs = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
#         output:
#             dict = os.path.join(OUTDIR,"{key}", "contigs.flt.dict")  
#         threads: threads_fn(rulename)
#         resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#         benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
#         log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}" + rulename
#         conda:
#             THIS_FILE_DIR / "envs/samtools.yaml"
#         shell:
#             "samtools dict {input} | cut -f1-3 > {output} 2> {log}"
#
#     # Generate bam files 
#     rule minimap:
#         input:
#             fw = read_fw,
#             rv = read_rv,
#             mmi = os.path.join(OUTDIR, "{key}", "contigs.flt.mmi"),
#             dict = os.path.join(OUTDIR,"{key}", "contigs.flt.dict"),
#             # fq = lambda wildcards: sample2path[wildcards.sample],
#             # mmi = os.path.join(OUTDIR,"contigs.flt.mmi"),
#             # dict = os.path.join(OUTDIR,"contigs.flt.dict")
#         output:
#             bam = bamfiles_before, 
#             # bam = temp(os.path.join(OUTDIR,"mapped/{sample}.bam"))
#         resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#         benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
#         log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
#         conda:
#             THIS_FILE_DIR / "envs/minimap2.yaml"
#         shell:
#             # See comment over rule "dict" to understand what happens here
#             "minimap2 -t {threads} -ax sr {input.mmi} {input.fw} {input.rv} -N 5"
#             " | grep -v '^@'"
#             " | cat {input.dict} - "
#             " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
#             " > {output.bam} 2> {log}"
#
# else:
#     # Run strobealign to get the abundances  
#     rulename = "Strobealign_bam_default"
#     rule Strobealign_bam_default:
#             input: 
#                 fw = read_fw,
#                 rv = read_rv,
#                 contig = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
#             output:
#                 bamfiles_before,
#             threads: threads_fn(rulename)
#             resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#             benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
#             log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
#             conda: THIS_FILE_DIR / "envs/strobe_env.yaml"
#             shell:
#                 """
#                 strobealign -t {threads} {input.contig} {input.fw} {input.rv} > {output} 2> {log}
#                 """
#
# Sort the bam files and index them
rulename="sort"
rule sort:
    input:
        # OUTDIR / "{key}/assembly_mapping_output/mapped/{id}.bam",
        bamfiles,
    output:
        OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    conda: THIS_FILE_DIR / "envs/samtools.yaml"
    shell:
        """
	samtools sort --threads {threads} {input} -o {output} 2> {log}
	# samtools index {output} 2>> {log}
	"""

rulename = "centrifuge_db"
rule centrifuge_db:
    input:
        table = "/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/centri/seqid2taxid.map",
        tax = "/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/centri/taxonomy/nodes.dmp",
        names = "/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/centri/taxonomy/names.dmp",
        contigs = "/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/centri/input-sequences.fna"
    output:
        "/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/centri/abv"
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    # default_target: True
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/centrifuge.yaml"
    shell:
        """
             centrifuge-build -p {threads} --conversion-table {input.table} \
            --taxonomy-tree {input.tax} --name-table {input.names} \
             {input.contigs} {output}
        """
rulename = "kraken_db"
rule kraken_db:
    threads: threads_fn(rulename)
    # default_target: True
    threads: 60
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/kraken2.yaml"
    shell:
        """
             kraken2-build --standard --db kraken_database_4 --use-ftp --threads {threads}
        """


## Include the specific rules for each tool
include: THIS_FILE_DIR / "snakemake_modules/vamb_default.smk"
include: THIS_FILE_DIR / "snakemake_modules/metabat.smk"
include: THIS_FILE_DIR / "snakemake_modules/comebin.smk"
include: THIS_FILE_DIR / "snakemake_modules/metadecoder.smk"
include: THIS_FILE_DIR / "snakemake_modules/semibin.smk"
include: THIS_FILE_DIR / "snakemake_modules/taxonomy_classifiers.smk"
include: THIS_FILE_DIR / "snakemake_modules/gunc.smk"
include: THIS_FILE_DIR / "snakemake_modules/va_extra_taxonomy_classifiers.smk"




rulename = "semibin_checkm"
rule semibin_checkm_2:
    input: 
        bin_dir = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/airways/PRJNA783873/before_may16_earlier_results/semibin/samples/output_prerecluster"
    output:
        outdir = directory(OUTDIR /  "soil_before_recluster_semibin"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fa' --database_path {params.database}
        """

import os

# Define the base directory for your data
BASE_DIR = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/data/PRJNA949880_aerosol_subset"

# Get a list of all sample IDs from your fastq files
# This assumes your fastq files are named like SRRXXXXXXX_1.fastq and SRRXXXXXXX_2.fastq
SAMPLES = sorted(list(set([os.path.basename(f).split('_')[0] for f in os.listdir(BASE_DIR) if f.endswith('.fastq')])))

# Define a function to get the forward and reverse reads for a given sample
# Modified to return a list of paths directly
def get_fastq_reads(wildcards):
    fw_read = os.path.join(BASE_DIR, f"{wildcards.sample}_1.fastq")
    rv_read = os.path.join(BASE_DIR, f"{wildcards.sample}_2.fastq")
    return [fw_read, rv_read] # Return a list of paths

rule all_spades:
    input:
        # This rule ensures that Snakemake knows what to build.
        # It will create an output directory for each sample processed by spades.
        # expand(os.path.join(OUTDIR, "spades_output", "{sample}"), sample=SAMPLES)
        bamfiles_before_out = expand(OUTDIR / "assembly_mapping_output/mapped/{sample}.bam", sample=SAMPLES)

# rulename = "spades"
# rule qc:
#     input:
#         get_fastq_reads
#     output:
#         read_fw_after_fastp = "data/sample_{sample}/reads_fastp/1.qc.fastq.gz",
#         read_rv_after_fastp =  "data/sample_{sample}/reads_fastp/2.qc.fastq.gz"
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     conda: THIS_FILE_DIR / "envs/fastp.yaml"
#     shell:
#             'fastp -i {input[0]} -I {input[1]} '
#             '-o {output.read_fw_after_fastp} -O {output.read_rv_after_fastp} --trim_tail1 1 --cut_tail --trim_poly_g --dedup --length_required 80 '
#             '--thread {threads} '

# rulename = "spades"
# rule spades_2:
#     input:
#         read_fw_after_fastp = "data/sample_{sample}/reads_fastp/1.qc.fastq.gz",
#         read_rv_after_fastp =  "data/sample_{sample}/reads_fastp/2.qc.fastq.gz"
#     output:
#         outdir = directory(os.path.join(OUTDIR, "spades_output", "{sample}"))
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     conda: THIS_FILE_DIR / "envs/spades_env.yaml"
#     shell:
#        "spades.py --meta "
#        "-t {threads} -k 21,33,55,77,99 "
#         "-o {output.outdir} -1 {input[0]} -2 {input[1]} "
#        "-t {threads} --memory {resources.mem_gb}" 


INDEX_SIZE = "12G" # get_config("index_size", "12G", r"[1-9]\d*[GM]$")
# Index resulting contig-file with minimap2
rulename = "minimap"

MIN_CONTIG_LEN = 2000

read_fw_after_fastp = "data/sample_{sample}/reads_fastp/1.qc.fastq.gz",
read_rv_after_fastp =  "data/sample_{sample}/reads_fastp/2.qc.fastq.gz"

# Cat the contigs together in one file to later map each pair of reads against all the contigs together
rule cat_contigs_2:
    input: 
        expand(os.path.join(OUTDIR, "spades_output", "{sample}", "contigs.fasta"), sample = SAMPLES),
    output: OUTDIR / "spades_output/contigs.flt.fna.gz"
    threads: threads_fn(rulename)
    params: script =  SRC_DIR / "concatenate.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell: 
        "python {params.script} {output} {input} -m {MIN_CONTIG_LEN} "  


rule index:
    input:
        contigs = OUTDIR / "spades_output/contigs.flt.fna.gz"    
    output:
        mmi = os.path.join(OUTDIR, "spades_output", "contigs.flt.mmi")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: 
        THIS_FILE_DIR / "envs/minimap2.yaml"
    shell:
        "minimap2 -I {INDEX_SIZE} -d {output} {input} "

# This rule creates a SAM header from a FASTA file.
# We need it because minimap2 for truly unknowable reasons will write
# SAM headers INTERSPERSED in the output SAM file, making it unparseable.
# To work around this mind-boggling bug, we remove all header lines from
# minimap2's SAM output by grepping, then re-add the header created in this
# rule.
rule dict:
    input:
        contigs = OUTDIR / "spades_output/contigs.flt.fna.gz"    
    output:
        dict = os.path.join(OUTDIR,"spades_output", "contigs.flt.dict")  
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda:
        THIS_FILE_DIR / "envs/minimap2.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} "

bamfiles_before = OUTDIR / "assembly_mapping_output/mapped/{sample}.bam"
# Generate bam files 
rule minimap:
    input:
        fw = read_fw_after_fastp,
        rv = read_rv_after_fastp,
        mmi = os.path.join(OUTDIR, "spades_output", "contigs.flt.mmi"),
        dict = os.path.join(OUTDIR,"spades_output", "contigs.flt.dict"),
    output:
        bam = bamfiles_before, 
        # bam = temp(os.path.join(OUTDIR,"mapped/{sample}.bam"))
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda:
        THIS_FILE_DIR / "envs/minimap2.yaml"
    shell:
        # See comment over rule "dict" to understand what happens here
        "minimap2 -t {threads} -ax sr {input.mmi} {input.fw} {input.rv} -N 5"
        " | grep -v '^@'"
        " | cat {input.dict} - "
        " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
        " > {output.bam} "
