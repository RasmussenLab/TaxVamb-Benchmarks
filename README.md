# Benchmark taxvamb against other binning tools

This is the code for the workflow for creating most of the benchmarks in the Taxvamb paper.

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

The workflows runs several taxonomy annotation tools and checkm, which all require databases.  
These databases should be installed and their paths set in the config file at config/config.yaml:

| Tool | Database Version | Description / Notes |
| :--- | :--- | :--- |
| **Metabuli** | GTDB v214.1 + T2T-CHM13v2.0 |Default database: Complete Genome/Chromosome, CheckM completeness > 90 and contamination < 5 + A human genome (T2T-CHM13v2.0). Installed using `metabuli databases` option|
| **MMseqs2** | GTDB v220 | Installed using the `mmseqs databases` option |
| **MMseqs2** | TrEMBL Release 2025_01 | Installed using the `mmseqs databases` option |
| **MMseqs2** | Kalmari (v3.7) | Installed using the `mmseqs databases` option |
| **Centrifuge** | NCBI RefSeq Release 229 | See https://www.ccb.jhu.edu/software/centrifuge/manual.shtml under 'Database download and index building' |
| **Kraken2** | Standard RefSeq (2024-12-28) | Pre-built index from [Langmead AWS](https://benlangmead.github.io/aws-indexes/k2). Includes archaea, bacteria, viral, plasmid, human, and UniVec_Core. |
| **Checkm** | Diamond db | Install using `checkm2 database --download`|

---

### Running the workflow

To dryrun the pipeline run
```
make benchmark_dryrun config=<bam_config_file>
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

##### Semibin
To use semibin with gpu in config/config.yaml set
`semibin_use_gpu: False` to True
Due to cuda problems we in this rule we use our HPC specific cuda modules: `module load cuda/12.2`. Therefore to get semibin to run on GPU you might need to play around with the 'semibinGPU' rule in snakemake_modules/semibin.smk


#### For longread datasets
For the longread benchmarks use semibin with the `--sequencing-type=long_read` flag.

---------------
### Running on a cluster with snakemake submiting jobs 
For running PlasMAAG on a cluster with snakemake submiting jobs see the documentation for snakemake [here](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html)  
An example is provided below for reference using slurm running PlasMAAG from reads:
Start off by installing the cluster-generic executor plugin for snakemake
```
pip install snakemake-executor-plugin-cluster-generic
```
Then run the PlasMAAG snakemake pipeline:
```
snakemake --use-conda --snakefile <path_to_snakefile> --config read_assembly_dir=<reads_and_assembly_dir_file> output_directory=<output_directory> \
  --jobs 2 --max-jobs-per-second 5 --max-status-checks-per-second 5 --latency-wait 60 \
  --executor cluster-generic --cluster-generic-submit-cmd 'sbatch --job-name {rule} --time={resources.walltime} --cpus-per-task {threads} --mem {resources.mem_gb}G'
```

#### Resources for the different snakemake rules when using snakemake directly
To define resources for the specific snakemake rules edit the `config/config.yaml` file
For more information see the ["Resources" section](#Resources).

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

