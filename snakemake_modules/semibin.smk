# rulename = "semibin"
# rule semibin:
#     input:
#         contigs = contigs_all,
#         bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
#     output:
#         semibin = directory(OUTDIR /  "{key}/semibin"),
#         semibin_bins = directory(OUTDIR /  "{key}/semibin/bins"),
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/semibin.yaml"
#     shell:
#         """
#         rm -rf {output.semibin}
#         SemiBin2 multi_easy_bin -i {input.contigs} -b {input.bamfiles} -o {output.semibin} \
#         --separator C -t {threads} --write-pre-reclustering-bins --self-supervised # &> {log}
#         """


rulename = "semibinGPU"
rule semibinGPU:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        semibin = directory(OUTDIR /  "{key}/semibin"),
        semibin_bins = directory(OUTDIR /  "{key}/semibin/bins"),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    # conda: THIS_FILE_DIR / "envs/semibin.yaml"
    params:
        env = "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/pixi"
    shell:
        """
        rm -rf {output.semibin}
        module load cuda/12.2
        CONDA_OVERRIDE_CUDA="12.2" pixi run --manifest-path {params.env} SemiBin2 \
        multi_easy_bin -i {input.contigs} -b {input.bamfiles} -o {output.semibin} \
        --separator C -t {threads} --engine gpu --write-pre-reclustering-bins --self-supervised # &> {log}
        """

rulename = "semibin"
rule semibin_LR:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        semibin = directory(OUTDIR /  "{key}/semibin_LR"),
        semibin_bins = directory(OUTDIR /  "{key}/semibin_LR/bins"),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/semibin.yaml"
    shell:
        """
        rm -rf {output.semibin}
        SemiBin2 multi_easy_bin --engine gpu --sequencing-type=long_read -i {input.contigs} -b {input.bamfiles} -o {output.semibin} \
        --separator C -t {threads} --write-pre-reclustering-bins --self-supervised # &> {log}
        """

rulename = "semibin_checkm"
rule semibin_checkm_LR:
    input: 
        bin_dir = OUTDIR /  "{key}/semibin_LR/bins",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/semibin_LR"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'gz' --database_path {params.database}
        """

rulename = "semibin_checkm"
rule semibin_checkm:
    input: 
        bin_dir = OUTDIR /  "{key}/semibin/bins",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/semibin"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'gz' --database_path {params.database}
        """
