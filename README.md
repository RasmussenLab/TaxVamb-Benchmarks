# Benchmark taxvamb against other binning tools

This is the reposityory for creating most of the benchmarks in the Taxvamb paper.

## Installation

Clone the repository and install the package using conda
```
git clone https://github.com/RasmussenLab/PlasMAAG
conda env create -n Benchmark_binners --file=taxvamb_paper_benchmarks/envs/benchmark_env.yaml
```
To use the program activate the conda environment
```
conda activate Benchmark_binners
```

Metabat is available as an docker image. To run it it is necessary to have singularity installed.
See [documentation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) for how to install 

### Running using snakemake CLI directly 
For using snakemake refer to the snakemake documentation: <https://snakemake.readthedocs.io/en/stable/>

#### Running from Reads using snakemake directly
To run the entire pipeline pass in a whitespace separated file using to the config flag in the snakemake CLI:

```
snakemake --use-conda --cores <number_of_cores> --snakefile <path_to_snakefile> --config bam_contig=<bam_contig_file> 
```

The <bam_contig_file> could look like:

``` 
sample 				bamfile				           contig
test                test_data/bam/sample_0.bam     test_data/contigs/contigs.fasta
test                test_data/bam/sample_1.bam     test_data/contigs/contigs.fasta
```
:heavy_exclamation_mark: Notice the header names are required to be: sample, bamfile and contig  
:heavy_exclamation_mark: Furthermore the contig file should just be the same for each row.

The input files are the following
1) A contig file, if using several samples this is all of them concatenated
2) BAM file(s) from mapping short reads to concatenated contig

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

#### Long read
For longread change the semibin commands to use it's LR mode.
