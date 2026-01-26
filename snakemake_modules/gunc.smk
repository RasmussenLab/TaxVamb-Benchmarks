

# name: [filepath, bin suffix]
all_bin_dirs = {
            # "reclustering_default_vamb": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/default_vamb", ".fna"],
            # Binners
            # "centrifuge_taxvamb_no_predictor": [OUTDIR / "{key}/formatted_vamb_bins_filtered/centrifuge_taxvamb_no_predictor", ".fna"],
            # "kraken_taxvamb_default": [OUTDIR / "{key}/formatted_vamb_bins_filtered/kraken_taxvamb_default", ".fna"],
            # "kraken_taxvamb_no_predictor": [OUTDIR / "{key}/formatted_vamb_bins_filtered/kraken_taxvamb_no_predictor", ".fna"],
            # "metabuli_taxvamb_default": [OUTDIR / "{key}/formatted_vamb_bins_filtered/metabuli_taxvamb_default", ".fna"],
            # "metabuli_taxvamb_no_predictor": [OUTDIR / "{key}/formatted_vamb_bins_filtered/metabuli_taxvamb_no_predictor", ".fna"],
            # "run_taxvamb_centrifuge": [OUTDIR / "{key}/formatted_vamb_bins_filtered/run_taxvamb_centrifuge", ".fna"],
            # "run_taxvamb_gtdb_no_predictor": [OUTDIR / "{key}/formatted_vamb_bins_filtered/run_taxvamb_gtdb_no_predictor", ".fna"],
            # "run_taxvamb_gtdb_w_unknown": [OUTDIR / "{key}/formatted_vamb_bins_filtered/run_taxvamb_gtdb_w_unknown", ".fna"],
            # "default_vamb": [OUTDIR / "{key}/vamb_default_bins_filtered", ".fna"],
            # "reclustering_centrifuge_taxvamb_no_predictor": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/centrifuge_taxvamb_no_predictor", ".fna"],
            # "reclustering_kraken_taxvamb_default": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/kraken_taxvamb_default", ".fna"],
            # "reclustering_kraken_taxvamb_no_predictor": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/kraken_taxvamb_no_predictor", ".fna"],
            # "reclustering_metabuli_taxvamb_default": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/metabuli_taxvamb_default", ".fna"],
            # "reclustering_metabuli_taxvamb_no_predictor": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/metabuli_taxvamb_no_predictor", ".fna"],
            # "reclustering_run_taxvamb_centrifuge": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/run_taxvamb_centrifuge", ".fna"],
            # "reclustering_run_taxvamb_gtdb_no_predictor": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/run_taxvamb_gtdb_no_predictor", ".fna"],
            # "reclustering_run_taxvamb_gtdb_w_unknown": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/run_taxvamb_gtdb_w_unknown", ".fna"],
            # "metabat": [OUTDIR / "{key}/metabat/metabat", ".fa"],
            # "metadecoder": [OUTDIR / "{key}/metadecoder/clusters", ".fasta"],
            "comebin": [OUTDIR / "{key}/comebin/comebin_res/comebin_res_bins", ".fa"],
            "semibin": [OUTDIR / "{key}/semibin/bins", ".fa"],
            # "kalmari_taxvamb_default": [OUTDIR / "{key}/formatted_vamb_bins_filtered/kalmari_taxvamb_default", ".fna"],
            # "trembl_taxvamb_default": [OUTDIR / "{key}/formatted_vamb_bins_filtered/trembl_taxvamb_default", ".fna"],
            # "trembl_taxvamb_default_reclustering": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/trembl_taxvamb_default", ".fna"],
            # "kalmari_taxvamb_default_reclustering": [OUTDIR / "{key}/reclustering/formatted_vamb_bins/kalmari_taxvamb_default", ".fna"],
            # # "semibin_LR": [OUTDIR / "{key}/semibin_LR/bins", ".fa"],
                }
bin_dir_names = all_bin_dirs.keys()

rulename = "gunc"
rule gunc:
    input:
        # semi_bin_decompressed = OUTDIR / "{key}/tmp/semibin_decompressed.done", # recomment me
        bin_dir = lambda wildcards: all_bin_dirs[wildcards.bin_dir][0]
    output: 
        dir = directory(OUTDIR / "{key}/gunc/{bin_dir}"),
        file = OUTDIR / "{key}/gunc/{bin_dir}/GUNC.progenomes_2.1.maxCSS_level.tsv"
    params:
        suffix = lambda wildcards: all_bin_dirs[wildcards.bin_dir][1],
        database = config.get("gunc_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/gunc.yaml"
    shell:
        """
        # Input dir should contain files in fasta format -- need to gzip for some
        rm -rf {output.dir}
        mkdir -p {output.dir}

        if [ "{wildcards.bin_dir}" == "semibin" ]; then
            gzip -dk {input.bin_dir}/*.gz 
        fi

        gunc run --input_dir {input.bin_dir} --file_suffix {params.suffix} --db_file {params.database} --threads {threads} --out_dir {output.dir} 
        """


rulename = "collect_gunc"
rule collect_gunc:
    input:
        expand(OUTDIR / "{{key}}/gunc/{bin_dir}",  bin_dir = bin_dir_names)
    output: 
        OUTDIR / "{key}/tmp/gunc.done"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/gunc.yaml"
    shell:
        """
        touch  {output}
        """
