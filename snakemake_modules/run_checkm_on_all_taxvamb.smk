## Collecting
all_bin_dirs_clas = {
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown/vaevae_clusters_split.tsv",
    "kraken_taxvamb_default": OUTDIR / "{key}/kraken_taxvamb_default/vaevae_clusters_split.tsv",
    "run_taxvamb_centrifuge": OUTDIR / "{key}/centrifuge_taxvamb/vaevae_clusters_split.tsv",
    "metabuli_taxvamb_default": OUTDIR / "{key}/metabuli_taxvamb_default/vaevae_clusters_split.tsv",
    "kalmari_taxvamb_default": OUTDIR / "{key}/kalmari_taxvamb_default/vaevae_clusters_split.tsv",
    "trembl_taxvamb_default": OUTDIR / "{key}/trembl_taxvamb_default/vaevae_clusters_split.tsv",
}
bin_dir_names_clas = all_bin_dirs_clas.keys()

# Prevent bins_clas wildcard from matching paths with slashes to avoid conflicts with reclustering outputs
wildcard_constraints:
    bins_clas="[^/]+"

rulename = "format_bins_class"
rule format_bins_class:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins = lambda wildcards: all_bin_dirs_clas[wildcards.bins_clas]
    output:
        directory = directory(OUTDIR / "{key}/formatted_vamb_bins_filtered/{bins_clas}"),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """


rulename = "checkm_class"
rule checkm_class:
    input: 
        bin_dir = directory(OUTDIR / "{key}/formatted_vamb_bins_filtered/{bins_clas}"),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/{bins_clas}"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

rulename = "collect_checkm"
rule collect_checkm:
    input:
        expand(OUTDIR /  "{{key}}/checkm2/{bins_clas}", bins_clas = bin_dir_names_clas)
    output: 
        OUTDIR / "{key}/tmp/checkm.done"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    shell:
        """
        touch  {output}
        """
