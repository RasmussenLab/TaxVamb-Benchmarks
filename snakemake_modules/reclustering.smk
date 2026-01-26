
rulename = "format_bins_class_recluster"
rule rename_vamb:
    input:
        bins = OUTDIR / "{key}/vamb_default/vae_clusters_split.tsv",
        latent = OUTDIR / "{key}/vamb_default/latent.npz",
    output:
        bins = OUTDIR / "{key}/vamb_default/vaevae_clusters_split.tsv",
        latent = OUTDIR / "{key}/vamb_default/vaevae_latent.npz",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    shell:
        """
            cp {input.bins} {output.bins}
            cp {input.latent} {output.latent}
        """

rulename = "recluster"
rule recluster:
    input: 
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = lambda wildcards: all_bin_dirs_recluster[wildcards.bins_recluster],
        bins = lambda wildcards: all_bin_dirs_recluster[wildcards.bins_recluster] / "vaevae_clusters_split.tsv",
        compo = lambda wildcards: all_bin_dirs_recluster[wildcards.bins_recluster] / "composition.npz",
    output:
        directory = directory(OUTDIR / "{key}/reclustering/{bins_recluster}/output"),
        headers = OUTDIR / "{key}/reclustering/{bins_recluster}/headers.txt",
        bins = OUTDIR / "{key}/reclustering/{bins_recluster}/output/clusters_reclustered.tsv",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    params:
        env = THIS_FILE_DIR / "reclustering"
    shell:
        """
        grep -E "^>" {input.contigs_decompressed} | cut -c 2- > {output.headers}
        rm -rf {output.directory}
        mkdir -p {output.directory}
        pixi run --manifest-path {params.env} start --num_process {threads} \
        {input.bins} {input.directory}/vaevae_latent.npz \
        {input.contigs_decompressed} \
        {output.headers} \
        {output.directory} \
        kmeans \
        """

# Run taxvamb 
rulename = "format_bins_class_recluster"
rule format_bins_class_recluster:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins = OUTDIR / "{key}/reclustering/{bins_recluster}/output/clusters_reclustered.tsv",
    output:
        directory = directory(OUTDIR / "{key}/reclustering/formatted_vamb_bins/{bins_recluster}"),
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


rulename = "checkm_class_recluster"
rule checkm_class_recluster:
    input:
        bin_dir = OUTDIR / "{key}/reclustering/formatted_vamb_bins/{bins_recluster}",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/reclustering/{bins_recluster}"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """
