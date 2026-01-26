# Benchmark taxvamb against other binning tools

This is the code for creating the benchmarks for the taxvamb paper for benchmarking the following binners:
- Taxvamb with different classifiers and databases:
  - Kalamari + GTDB
  - MMseqs + GTDB / TrEMBL/ Kalmari
  - Centrifuge + NCBI RefSEQ
  - Kraken2 + NCBI RefSEQ
- Metabat
- Semibin2
- Comebin
- Metadecoder
Additionally, each binner is assesed using checkm2 and GUNC - with all vamb derrived binners being assesed before and after reclustering.  

The orginal workflow was ran on the ESRUM cluster: https://cbmr-data.github.io/esrum/overview.html


## Installation

Clone the repository and install the package using conda

```
git clone https://github.com/las02/taxvamb_paper_benchmarks
conda env create -n Benchmark_binners --file=taxvamb_paper_benchmarks/envs/benchmark_env.yaml
```
To use the program activate the conda environment
```
conda activate Benchmark_binners
```

Metabat is available as an docker image. To run it it is necessary to have singularity installed.
See [documentation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) for how to install the software

Taxconverter should be installed to a conda environemnt, see documentation for installation [taxconverter](https://github.com/RasmussenLab/taxconverter/tree/fix_taxonomy_assignment_for_unclassified_vs_no_rank)


The workflows runs several taxonomy annotation tools and checkm, which all require databases.  
These databases should be installed and their paths set in the config file at `config/config.yaml`:

| Tool | Database Version | Description / Notes |
| :--- | :--- | :--- |
| **Metabuli** | GTDB v214.1 + T2T-CHM13v2.0 |Default database: Complete Genome/Chromosome, CheckM completeness > 90 and contamination < 5 + A human genome (T2T-CHM13v2.0). Installed using `metabuli databases` option|
| **MMseqs2** | GTDB v220 | Installed using the `mmseqs databases` option |
| **MMseqs2** | TrEMBL Release 2025_01 | Installed using the `mmseqs databases` option |
| **MMseqs2** | Kalmari (v3.7) | Installed using the `mmseqs databases` option |
| **Centrifuge** | NCBI RefSeq Release 229 | See https://www.ccb.jhu.edu/software/centrifuge/manual.shtml under 'Database download and index building' |
| **Kraken2** | RefSeq (2024-12-28) | Pre-built index from [Langmead AWS](https://benlangmead.github.io/aws-indexes/k2). Includes archaea, bacteria, viral, plasmid, human, and UniVec_Core. |
| **Checkm** | Diamond db | Install using `checkm2 database --download`|

---

Lastly Taxconverter should then be installed in a conda environment called taxconv, see the documentation for installation [Taxconverter documentation](https://github.com/RasmussenLab/taxconverter/tree/fix_taxonomy_assignment_for_unclassified_vs_no_rank). 

### Running the workflow

#### Running the entire workflow

To dry the workflow there is several options
```
make benchmark_dryrun config=<bam_config_file>     # Dry-running the workflow

make benchmark_run config=<bam_config_file>        # Running the workflow locally

make benchmark_run_slurm config=<bam_config_file>  # Running the workflow on SLURM
```

The <bam_contig_file> should look like:

``` 
sample 				bamfile				           contig
test                test_data/bam/sample_0.bam     test_data/contigs/contigs.fasta
test                test_data/bam/sample_1.bam     test_data/contigs/contigs.fasta
```
:heavy_exclamation_mark: The header names are required to be: sample, bamfile and contig  
:heavy_exclamation_mark: The contig file should be the same for each row.

The input files are the following
1) A contig file of all contigs concatenated. If using several samples this is all of them concatenated.
2) BAM file(s) from mapping short reads to the concatenated contig file

To generate the BAM file(s) and the contig files from assemblies and reads use the `map_snakefile.smk` pipeline.
The snakefile can be run using:

```
make benchmark_dryrun config=<read_assembly_dir>     # Dry-running the workflow

make benchmark_run_slurm config=<read_assembly_dir>  # Running the workflow on SLURM
```

Where read_assembly_dir should have the following structure:

```
sample	   read1                          read2                         contig                                           
sample_1   im/a/path/to/sample_1/read1    im/a/path/to/sample_1/read2   path/sample_1/contig.fasta
sample_2   im/a/path/to/sample_2/read1    im/a/path/to/sample_2/read2   path/sample_2/contig.fasta          
```

#### Running specific tools only

For running a specific tool snakemake can be called alone with a target output file eg.  

```
  snakemake --snakefile snakefile.smk -c 100 --software-deployment-method apptainer --use-conda  --config bam_contig=<bam_contig_file> <target_output_files> 
```
For more information see the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/).  
To figure out which files, need to be written refer to the target rule (all) in the `snakemake.smk`.  
For now it looks like:  

```
rule all:
    input:
        checkm_semibin = expand(OUTDIR /  "{key}/checkm2/semibin", key=sample_id.keys()),
        checkm_comebin = expand(OUTDIR /  "{key}/checkm2/comebin", key=sample_id.keys()),
        checkm_metadecoder = expand(OUTDIR /  "{key}/checkm2/metadecoder", key=sample_id.keys()),
        checkm_metabat = expand(OUTDIR /  "{key}/checkm2/metabat", key=sample_id.keys()),
        checkm_default_vamb = expand(OUTDIR /  "{key}/checkm2/default_vamb",key=sample_id.keys()),
        checkm2_taxvamb = expand(OUTDIR /  "{key}/tmp/checkm.done",key=sample_id.keys()), 
        gunc = expand(OUTDIR / "{key}/tmp/gunc.done", key=sample_id.keys()),
```
So to run semibin and checkm for a <bam_contig_file> which look like:  
``` 
sample 				bamfile				           contig
test                test_data/bam/sample_0.bam     test_data/contigs/contigs.fasta
test                test_data/bam/sample_1.bam     test_data/contigs/contigs.fasta
```
you would set the target file as `test/checkm2/semibin`  

To only run gunc for some tools, it is a few more steps. Here you need to  
- 1) set the target rule (for the above example) as `test/tmp/gunc.done`
- 2) in the `snakemake_modules/gunc.smk` files `all_bin_dirs` variable comment out the output files you are not interested in.  
For example for only running comebin and semibin this variable would look like:
```
all_bin_dirs = {
            "comebin": [OUTDIR / "{key}/comebin/comebin_res/comebin_res_bins", ".fa"],
            "semibin": [OUTDIR / "{key}/semibin/bins", ".fa"],
	    ##  Rest commented out ...
	    ##  ...
                }
```
The same applies if you only want to run some of the taxvamb version (ie. with a specific database), here you should have the following target file `test/tmp/checkm.done` and in `snakemake_modules/run_checkm_on_all_taxvamb.smk`, select the target files, ie. in the example only run taxvamb with the gtdb database:
```
all_bin_dirs_clas = {
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown/vaevae_clusters_split.tsv",
    ##  Rest commented out ...
    ##  ...
}
```

### Resources 
The pipeline can be configurated in: ``` config/config.yaml ```

Here, the resources for each rule can be configurated as follows
```
spades:
  walltime: "15-00:00:00"
  threads: 16
  mem_gb: 60
```
if no resources are configurated for a rule the defaults will be used which are also defined in: ``` config/config.yaml ```  as
```
default_walltime: "48:00:00"
default_threads: 16
default_mem_gb: 50
```
If these exceed the resources available they will be scaled down to match the hardware available. 

#### Using GPU
##### Taxvamb and Vamb
To use taxvamb and vamb with GPU assign the job to a GPU node. This can for example be done in ``` config/config.yaml ``` for the rules by adding `gpu: " --partition=gpuqueue --gres=gpu:1 "` on the HPC used for orginally running the workflow:
```
run_taxvamb_kraken:
  walltime: "20-00:00:00"
  mem_gb: 500
  threads: 64
  gpu: " --partition=gpuqueue --gres=gpu:1 "
```
Additionally set `vamb_use_gpu: True` in the conifg file

##### Semibin
To use semibin with gpu in config/config.yaml set
`semibin_use_gpu: False` to True
Due to cuda problems we in this rule we use our HPC specific cuda modules: `module load cuda/12.2`. Therefore to get semibin to run on GPU you might need to play around with the 'semibinGPU' rule in snakemake_modules/semibin.smk. 

#### For longread datasets
For the longread benchmarks use semibin with the `--sequencing-type=long_read` flag.  
Additionally, for reclustering the dbscan algorithm shoud be used instead of kmeans.
Lastly, for minimap `-ax map-hifi` should be used.

#### Running taxvamb/vamb without predictor
For the taxvamb/vamb runs pass in the `--no_predictor` flag

### Tools which crashed internally, and alternative ways of running them - along with the resourses used for running the tools.
For benchmarking for figure 3 in the paper the following 4 runs crashed internally.

#### Semibin
Semibin(v.2.1.0) crashes in the Vaginal and the Salvia samples.
- For the Vaginal sample Semibin does not find any bins in one of the samples which causes downstream steps to crash (`/log_files_for_crashed_runs/Vaginal_SemiBin_v2.1.0.log` ). Semibin version 2.2.0 should fix this issue (https://github.com/BigDataBiology/SemiBin/releases/tag/v2.2.0), but does not change performance of the tool (according to patchnotes). We therefore ran semibin version 2.2.0 on this dataset, but we still got the same error (`/log_files_for_crashed_runs/Vaginal_SemiBin_v2.2.0.log`). This is due to another bug with no bins found using the `--write-pre-reclustering-bins` flag , removing this flag (as we don't need the pre-reclustering bins for this sample) and upgrading to the newest version of semibin fixes the issue. 
- For the the Salvia dataset we get the following error described in: https://github.com/BigDataBiology/SemiBin/issues/211 and https://github.com/BigDataBiology/SemiBin/issues/201. There does not seem to be a fix for the issue, although the maintainer seems to be looking into it. See `/log_files_for_crashed_runs/Salvia_SemiBin.log` for logfiles

#### Comebin
Comebin (v1.0.3) crashes in the Human Gut (IBS) and Forest Soil samples.
In both datasets the error is described in the following Github issue: https://github.com/ziyewang/COMEBin/issues/17,
The logfiles for these runs can be found in: `/log_files_for_crashed_runs/Human_gut_IBS_ComeBin.log` and `/log_files_for_crashed_runs/Forest_soil_Comebin.log`
Here we instead ran comebin in single-sample mode. This is equivalent to in the pipeline in the `sample` column of the config files, assigning each read pair and their corresponding contig file to a different sample name.

---

## Taxonomy Formatting Scripts

The pipeline includes two Python scripts that convert taxonomy classifier outputs to a standardized format compatible with Taxvamb:

### format_metabuli.py

**Purpose**: Converts Metabuli taxonomy classification output to Taxvamb-compatible format.

**What it does**:
- Parses the Metabuli report file to build a hierarchical taxonomy structure mapping tax IDs to full lineages
- Converts the Metabuli classification file (contig → taxID) to use full taxonomy lineage strings
- Outputs tab-separated format: `contigs\tpredictions`

**Pipeline usage**: Called in Metabuli-based Taxvamb runs (see `snakemake_modules/taxvamb_using_mmseqs_classifications.smk`)

**Example**:
```bash
python format_metabuli.py <classification_file> <report_file> > formatted_taxonomy.tsv
```

### format_trembl_kalmari.py

**Purpose**: Standardizes and validates MMseqs2 taxonomy output from TrEMBL or Kalmari databases for Taxvamb.

**What it does**:
- **Fixes missing taxonomic levels**: Ensures all 7 standard ranks (k,p,c,o,f,g,s) are present by filling gaps with placeholder names derived from the nearest known parent taxon (e.g., `LEVEL_4_ADDED_FROM_g_Escherichia` for a missing family)
- **Filters non-standard entries**: Removes intermediate taxonomy entries prefixed with `-_` while preserving domain-level classifications and subspecies annotations
- **Validates output**: Confirms the formatted taxonomy has the correct number of levels for the deepest rank present

**Pipeline usage**: Called in two rules in `snakemake_modules/taxvamb_using_mmseqs_classifications.smk`:
- `run_taxvamb_kalmari`: Processes MMseqs2 output from Kalmari database (column 5)
- `run_taxvamb_trembl`: Processes MMseqs2 output from TrEMBL database (column 9)

**Example**:
```bash
cut -f1,5 mmseqs_output.tsv > cut.tsv
echo -e "contigs\tpredictions" > formatted.tsv
python format_trembl_kalmari.py cut.tsv >> formatted.tsv
```

**Why these scripts are needed**: Different taxonomy classifiers (Metabuli, MMseqs2 with TrEMBL/Kalmari) produce outputs in varying formats with inconsistent taxonomic hierarchies. These scripts standardize the outputs to ensure Taxvamb receives taxonomy data in a consistent, validated format regardless of the upstream classifier used.
The utilities for these conversions for the metabuli annotations can be found in the taxconverter tool, although the taxconverter tool does not support the conversions from teh trembl and kalmari database as these are not the recommened databases - but are added in this pipeline to assess the difference annotations make wrt. to the performance of taxvamb.

