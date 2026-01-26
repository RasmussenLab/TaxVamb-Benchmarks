rulename = "comebin"
rule comebin:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output: 
        comebin = directory(OUTDIR /  "{key}/comebin"),
        bins = directory(OUTDIR / "{key}/comebin/comebin_res/comebin_res_bins"),
    params:
        bamfiles_dir = str(OUTDIR) + "/{key}/assembly_mapping_output/mapped_sorted",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename), load=1
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/comebin.yaml"
    #conda: "GPUcomebin"
    shell:
        """
	rm -rf {output.comebin}
        run_comebin.sh -a {input.contigs_decompressed} \
        -o {output.comebin} \
        -p {params.bamfiles_dir} \
        -t {threads} # &> {log}
        """

rulename = "comebin_checkm"
rule comebin_checkm:
    input: 
        bin_dir = OUTDIR /  "{key}/comebin/comebin_res/comebin_res_bins",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/comebin"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fa' --database_path {params.database}
        """
