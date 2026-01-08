benchmark_dryrun: 
	snakemake -n -p --keep-going --snakefile snakefile.smk -c 100 --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --config bam_contig=$(config) --keep-incomplete all
benchmark_run:
	snakemake -p --keep-going --snakefile snakefile.smk -c 100 --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --config bam_contig=$(config) --keep-incomplete  all
benchmark_run_slurm:
	snakemake -p --keep-going  --snakefile snakefile.smk --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --keep-incomplete --latency-wait 60 --config bam_contig=$(config) --jobs 35 --max-jobs-per-second 35 --max-status-checks-per-second 35 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --job-name {rule}  --time={resources.walltime} --cpus-per-task {threads} --mem {resources.mem_gb}G {resources.gpu}"  all
map_dryrun:
	snakemake -n -p --keep-going --snakefile map_snakefile.smk -c 100 --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --config read_assembly_dir=$(config) --keep-incomplete all
map_runn:
	snakemake -p --keep-going --snakefile map_snakefile.smk -c 100 --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --config read_assembly_dir=$(config) --keep-incomplete  all
map_run_slurm:
	snakemake -p --keep-going  --snakefile map_snakefile.smk --nolock --rerun-triggers mtime --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --latency-wait 60 --keep-incomplete --config read_assembly_dir=$(config) --jobs 10 --max-jobs-per-second 10 --max-status-checks-per-second 10 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --job-name {rule}  --time={resources.walltime} --cpus-per-task {threads} --mem {resources.mem_gb}G {resources.gpu}"  all

