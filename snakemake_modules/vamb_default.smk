# Run vamb
rulename = "run_vamb"
rule run_vamb:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        directory = directory(OUTDIR / "{key}/vamb_default"),
        bins = OUTDIR / "{key}/vamb_default/vae_clusters_split.tsv",
        compo = OUTDIR / "{key}/vamb_default/composition.npz",
        latent = OUTDIR / "{key}/vamb_default/latent.npz",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory} # clean up dir eg. for failed runs
        vamb bin default {vamb_extra_arg} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} &> {log}
        """

# Format bins as fasta files as checkm require this. Filter out bins smaller than 200_000 bp
rulename = "format_bins_default_vamb_filtered"
rule format_bins_default_vamb_filtered:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins = OUTDIR / "{key}/vamb_default/vae_clusters_split.tsv",
    output:
        directory = directory(OUTDIR / "{key}/vamb_default_bins_filtered"),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory}
        """

rulename = "default_vamb_checkm"
rule default_vamb_checkm:
    input:
        bin_dir = OUTDIR / "{key}/vamb_default_bins_filtered",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/default_vamb"),
        file = OUTDIR /  "{key}/checkm2/default_vamb/quality_report.tsv",
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """



