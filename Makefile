benchmark_dryrun:
	snakemake -n -p --keep-going --snakefile snakefile.smk -c 100 --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --keep-incomplete --config bam_contig=$(config) output_directory=$(output) --directory $(output) 
benchmark_run:
	snakemake -p --keep-going --snakefile snakefile.smk -c 100 --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --keep-incomplete --config bam_contig=$(config) output_directory=$(output) --directory $(output) 
benchmark_run_slurm:
	snakemake -p --keep-going  --snakefile snakefile.smk --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --keep-incomplete --config bam_contig=$(config) output_directory=$(output) --directory $(output) --jobs 5 --max-jobs-per-second 5 --max-status-checks-per-second 5 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --job-name {rule}  --time={resources.walltime} --cpus-per-task {threads} --mem {resources.mem_gb}G {resources.gpu}" 
