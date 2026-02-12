# Taxvamb Benchmarking Pipeline

A Snakemake workflow for benchmarking metagenomic binning tools, developed for the Taxvamb paper.

## Overview

This pipeline benchmarks the following binning tools:

**Taxvamb** (with multiple taxonomy classifiers and databases):
- Metabuli + GTDB
- MMseqs2 + GTDB / TrEMBL / Kalamari
- Centrifuge + NCBI RefSeq
- Kraken2 + NCBI RefSeq

**Other binners:**
- Vamb (default)
- MetaBAT2
- SemiBin2
- COMEBin
- MetaDecoder

All binners are assessed using CheckM2 and GUNC. Vamb and TaxVamb are evaluated both before and after applying reclustering.

> **Note:** This workflow was originally developed and tested on the [ESRUM cluster](https://cbmr-data.github.io/esrum/overview.html).


## Installation

### Prerequisites

1. **Conda/Mamba** - For environment management
2. **Singularity/Apptainer** - Required for MetaBAT2 (runs as Docker image). See the [Singularity installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

### Setup

Clone the repository and create the conda environment:

```bash
git clone https://github.com/RasmussenLab/TaxVamb-Benchmarks
conda env create -n Benchmark_binners --file=TaxVamb-Benchmarks/envs/benchmark_env.yaml
conda activate Benchmark_binners
```

### Additional Dependencies

For running some of the TaxVamb Benchmarks **Taxconverter** should also be installed in the conda environment see [Taxconverter documentation](https://github.com/RasmussenLab/taxconverter/tree/a32d4fc) for installation instructions (link is to the commit used in the pipeline).


### Database Requirements

The pipeline requires several databases. Install them and configure their paths in `config/config.yaml` (Notice: install all databases takes around ~1TB of storage):

| Tool           | Database Version            | Description / Notes |
|:---------------|:----------------------------| :--- |
| **Metabuli**  | GTDB v214.1 + T2T-CHM13v2.0 | Default database: Complete Genome/Chromosome, CheckM completeness > 90 and contamination < 5 + human genome (T2T-CHM13v2.0). Install using `metabuli databases` |
| **MMseqs2**   | GTDB v220                   | Install using `mmseqs databases` |
| **MMseqs2**    | TrEMBL Release 2025_01      | Install using `mmseqs databases` |
| **MMseqs2**    | Kalamari (v3.7)             | Install using `mmseqs databases` |
| **Centrifuge** | NCBI RefSeq Release 229     | See [Centrifuge manual](https://www.ccb.jhu.edu/software/centrifuge/manual.shtml) under "Database download and index building" |
| **Kraken2**    | RefSeq (2024-12-28)         | Pre-built index from [Langmead AWS](https://benlangmead.github.io/aws-indexes/k2). Includes archaea, bacteria, viral, plasmid, human, and UniVec_Core |
| **CheckM2**    | Diamond db                  | Install using `checkm2 database --download` |
| **GUNC**       | gunc db                     | Install using `gunc download_db` |
--- 

### Common problems

On some systems when running the pipeline you might get:
```
/usr/bin/bash: conda: command not found
WorkflowError:
Error running conda info. Is conda installed and accessible? Error: Command 'conda info --json' returned non-zero exit status 127.
```
To fix this install conda inside the conda environment:
```
conda install conda
```

## Running the Workflow

### Quick Start

```bash
# Dry-run (preview what will be executed)
make benchmark_dryrun config=data_configs/example.tsv

# Run locally
make benchmark_run config=data_configs/example.tsv

# Run on SLURM cluster
make benchmark_run_slurm config=data_configs/example.tsv
```

### Input Configuration File

The configuration file is a tab-separated file with three required columns:

| Column | Description |
| :--- | :--- |
| `sample` | Sample identifier (same for all BAM files from one sample) |
| `bamfile` | Path to BAM file from read mapping |
| `contig` | Path to concatenated contig file |

**Example (`data_configs/example.tsv`):**

```tsv
sample              bamfile                        contig
test                test_data/bam/sample_0.bam     test_data/contigs/contigs.fasta
test                test_data/bam/sample_1.bam     test_data/contigs/contigs.fasta
```

> ⚠️ **Important:**
> - Header names must be exactly: `sample`, `bamfile`, and `contig`
> - The contig file path must be identical for all rows belonging to the same sample

1. **Contig file**: A FASTA file containing all contigs. For multi-sample datasets, concatenate all sample assemblies into a single file.
2. **BAM file(s)**: Alignment files from mapping short reads to the concatenated contig file.

To generate the BAM file(s) and the contig files from assemblies and reads see next section:

### Generating BAM Files from Raw Reads and Assemblie

If you have raw reads and assemblies, use the mapping pipeline to generate BAM files:

```bash
# Dry-run
make map_dryrun config=<read_assembly_config>

# Run on SLURM
make map_run_slurm config=<read_assembly_config>
```

The read/assembly configuration file format is the following:

```tsv
sample      read1                              read2                              contig
sample_1    path/to/sample_1/read1.fastq.gz    path/to/sample_1/read2.fastq.gz    path/to/sample_1/assembly.fasta
sample_2    path/to/sample_2/read1.fastq.gz    path/to/sample_2/read2.fastq.gz    path/to/sample_2/assembly.fasta
```

### Running Specific Tools

To run only specific tools, invoke Snakemake directly with target output files:

```bash
snakemake --snakefile snakefile.smk -c 100 \
    --software-deployment-method apptainer --use-conda \
    --config bam_contig=<config_file> <target_output_files>
```

For detailed Snakemake options, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

#### Available Target Rules

The main target rule (`all`) in `snakefile.smk` defines available outputs:

```python
rule all:
    input:
        checkm_semibin = expand(OUTDIR / "{key}/checkm2/semibin", key=sample_id.keys()),
        checkm_comebin = expand(OUTDIR / "{key}/checkm2/comebin", key=sample_id.keys()),
        checkm_metadecoder = expand(OUTDIR / "{key}/checkm2/metadecoder", key=sample_id.keys()),
        checkm_metabat = expand(OUTDIR / "{key}/checkm2/metabat", key=sample_id.keys()),
        checkm_default_vamb = expand(OUTDIR / "{key}/checkm2/default_vamb", key=sample_id.keys()),
        checkm2_taxvamb = expand(OUTDIR / "{key}/tmp/checkm.done", key=sample_id.keys()),
        gunc = expand(OUTDIR / "{key}/tmp/gunc.done", key=sample_id.keys()),
```

#### Example: Running SemiBin + CheckM2 Only

For a config file with `sample=test`, set the target file as:

```bash
snakemake ... test/checkm2/semibin
```

#### Running GUNC for Specific Tools

1. Set the target rule to `{sample}/tmp/gunc.done`
2. Edit `snakemake_modules/gunc.smk` and modify the `all_bin_dirs` variable to include only desired tools:

```python
all_bin_dirs = {
    "comebin": [OUTDIR / "{key}/comebin/comebin_res/comebin_res_bins", ".fa"],
    "semibin": [OUTDIR / "{key}/semibin/bins", ".fa"],
    # Comment out other tools...
}
```

#### Running Specific Taxvamb Configurations

1. Set the target file to `{sample}/tmp/checkm.done`
2. Edit `snakemake_modules/run_checkm_on_all_taxvamb.smk` and modify `all_bin_dirs_clas`:

```python
all_bin_dirs_clas = {
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown/vaevae_clusters_split.tsv",
    # Comment out other configurations...
}
```

#### Running Reclustering for Specific Binners

Use the `all_reclustering` target and modify `all_bin_dirs_recluster` in `snakefile.smk`:

```bash
snakemake --snakefile snakefile.smk -c 100 \
    --config bam_contig=<config_file> all_reclustering
```

```python
all_bin_dirs_recluster = {
    "default_vamb": OUTDIR / "{key}/vamb_default",
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown",
    # Comment out other binners...
}
```

---

## Configuration

### Resource Configuration

All resources are configured in `config/config.yaml`. Per-rule resources override defaults:

```yaml
# Default resources (used when rule-specific values are not defined)
default_walltime: "48:00:00"
default_threads: 16
default_mem_gb: 50

# Rule-specific resources
spades:
  walltime: "15-00:00:00"
  threads: 16
  mem_gb: 60
```

> **Note:** If specified resources exceed available hardware, they are automatically scaled down.

### GPU Configuration

#### Taxvamb and Vamb

1. Set `vamb_use_gpu: True` in `config/config.yaml`
2. Add GPU partition settings to the relevant rules:

```yaml
run_taxvamb_kraken:
  walltime: "20-00:00:00"
  mem_gb: 500
  threads: 64
  gpu: " --partition=gpuqueue --gres=gpu:1 "
```

#### SemiBin

1. Set `semibin_use_gpu: True` in `config/config.yaml`
2. The `semibinGPU` rule uses HPC-specific CUDA modules (`module load cuda/12.2`). You may need to adjust this in `snakemake_modules/semibin.smk` for your cluster.

### Long-Read Datasets

For long-read data:
- Use SemiBin with `--sequencing-type=long_read` flag
- Use DBSCAN algorithm instead of k-means for reclustering
- Use minimap2 with `-ax map-hifi` for read mapping

### Running Without Predictor

For Taxvamb/Vamb runs without the predictor, add the `--no_predictor` flag.

---

## Crashed runs

Some tools (SemiBin, COMEBin, GUNC) crashed internally on certain datasets. See `/log_files_for_crashed_runs/README.md` for details on root causes, workarounds, resources used, and log files.

---

## Taxonomy Formatting Scripts

The pipeline includes Python scripts to convert taxonomy classifier outputs to a standardized Taxvamb-compatible format.

### format_metabuli.py

Converts Metabuli classification output to Taxvamb format.

**Functionality:**
- Parses the Metabuli report file to build a hierarchical taxonomy structure (taxID → full lineage)
- Converts classification file (contig → taxID) to full taxonomy lineage strings
- Outputs tab-separated format: `contigs\tpredictions`

**Usage:**
```bash
python format_metabuli.py <classification_file> <report_file> > formatted_taxonomy.tsv
```

**Pipeline integration:** Called in Metabuli-based Taxvamb runs (see `snakemake_modules/taxvamb_using_mmseqs_classifications.smk`)

### format_trembl_kalmari.py

Standardizes MMseqs2 taxonomy output from TrEMBL or Kalamari databases.

**Functionality:**
- **Fixes missing taxonomic levels:** Ensures all 7 standard ranks (k, p, c, o, f, g, s) are present by filling gaps with placeholder names (e.g., `LEVEL_4_ADDED_FROM_g_Escherichia`)
- **Filters non-standard entries:** Removes intermediate taxonomy entries prefixed with `-_` while preserving domain-level classifications and subspecies annotations
- **Validates output:** Confirms correct number of taxonomy levels

**Usage:**
```bash
cut -f1,5 mmseqs_output.tsv > cut.tsv
echo -e "contigs\tpredictions" > formatted.tsv
python format_trembl_kalmari.py cut.tsv >> formatted.tsv
```

**Pipeline integration:** Used in `snakemake_modules/taxvamb_using_mmseqs_classifications.smk`:
- `run_taxvamb_kalmari`: Processes MMseqs2 Kalamari output (column 5)
- `run_taxvamb_trembl`: Processes MMseqs2 TrEMBL output (column 9)

### Why These Scripts Are Needed

Different taxonomy classifiers produce outputs in varying formats with inconsistent hierarchical structures. These scripts standardize outputs to ensure Taxvamb receives consistent, validated taxonomy data regardless of the upstream classifier.

> **Note:** The [Taxconverter](https://github.com/RasmussenLab/taxconverter) tool handles Metabuli conversions, but does not support TrEMBL or Kalamari databases. These are included in this pipeline to assess how different annotations affect Taxvamb performance.

