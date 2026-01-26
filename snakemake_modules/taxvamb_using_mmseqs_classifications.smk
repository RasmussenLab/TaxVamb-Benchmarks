rulename = "mmseqs2_kalmari"
rule mmseqs2_kalmari:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        mmseqs2 = directory(OUTDIR /  "{key}/classifiers/mmseqs2/kalmari"),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/kalmari_lca.tsv",
        tmp = directory(OUTDIR /  "{key}/classifiers/mmseqs2/tmp_kalmari"),
    threads: threads_fn(rulename)
    params:
        db = config.get("mmseqs_db_kalmari"),
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/mmseqs2.yaml"
    shell:
        """
        rm -rf {output.mmseqs2}
	mmseqs easy-taxonomy {input.contigs_decompressed} {params.db} {output.mmseqs2} {output.tmp} --search-type 3 --tax-lineage 1 
        mkdir -p {output.mmseqs2}
        """

rulename = "mmseqs2_gtdb"
rule mmseqs2_gtdb:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        mmseqs2 = directory(OUTDIR /  "{key}/classifiers/mmseqs2/gtdb"),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/gtdb_lca.tsv",
        tmp = directory(OUTDIR /  "{key}/classifiers/mmseqs2/tmp_gtdb"),
    threads: threads_fn(rulename)
    params:
        db = config.get("mmseqs_db_gtdb"),
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/mmseqs2.yaml"
    shell:
        """
        rm -rf {output.mmseqs2}
	mmseqs easy-taxonomy {input.contigs_decompressed} {params.db} {output.mmseqs2} {output.tmp} --tax-lineage 1 
        mkdir -p {output.mmseqs2}
        """


rulename = "mmseqs2_trembl"
rule mmseqs2_trembl:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        mmseqs2 = directory(OUTDIR /  "{key}/classifiers/mmseqs2/trembl"),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/trembl_lca.tsv",
        tmp = directory(OUTDIR /  "{key}/classifiers/mmseqs2/tmp_trembl"),
    threads: threads_fn(rulename)
    params:
        db = config.get("mmseqs_db_trembl"),
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/mmseqs2.yaml"
    shell:
        """
        rm -rf {output.mmseqs2}
	mmseqs easy-taxonomy {input.contigs_decompressed} {params.db} {output.mmseqs2} {output.tmp} --tax-lineage 1 
        mkdir -p {output.mmseqs2}
        """

# Run taxvamb 
rulename = "run_taxvamb_kalmari"
rule run_taxvamb_kalmari:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/kalmari_lca.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kalmari_taxvamb_default')),
        bins = OUTDIR / "{key}/kalmari_taxvamb_default/vaevae_clusters_split.tsv",
        compo = os.path.join(OUTDIR, '{key}','kalmari_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    params: script =THIS_FILE_DIR / "files_used_in_snakemake_workflow/format_trembl_kalmari.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        pip install typer

        # Make sure the taxonomy files produced by kalmari are in the format taxvamb expects
        cut -f1,5 {input.mmseqs2_out} > {input.mmseqs2_out}.cut.tsv
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        python {params.script} {input.mmseqs2_out}.cut.tsv >> {input.mmseqs2_out}.formatted

        rm -rf {output.directory}
        vamb bin taxvamb {vamb_extra_arg} --taxonomy {input.mmseqs2_out}.formatted --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """


# Run taxvamb 
rulename = "run_taxvamb_gtdb"
rule run_taxvamb_gtdb_w_unknown:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/gtdb_lca.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default_w_unknown')),
        bins = OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown/vaevae_clusters_split.tsv",
        compo = os.path.join(OUTDIR, '{key}','gtdb_taxvamb_default_w_unknown/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        # Format taxonomy file
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted_w_unknown
        sed 's/\t$/\tunknown/' {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted_w_unknown  # Taxonomy levels not classified should output 'unknown' instaed of an empty string ('')
        rm -rf {output.directory}
        vamb bin taxvamb {vamb_extra_arg} --taxonomy {input.mmseqs2_out}.formatted_w_unknown --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb 
rulename = "run_taxvamb_trembl"
rule run_taxvamb_trembl:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/trembl_lca.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'trembl_taxvamb_default')),
        bins = OUTDIR / "{key}/trembl_taxvamb_default/vaevae_clusters_split.tsv",
        compo = os.path.join(OUTDIR, '{key}','trembl_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    params: script =THIS_FILE_DIR / "files_used_in_snakemake_workflow/format_trembl_kalmari.py"
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        cut -f1,9 {input.mmseqs2_out} > {input.mmseqs2_out}.cut.tsv
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        python {params.script} {input.mmseqs2_out}.cut.tsv >> {input.mmseqs2_out}.formatted

        rm -rf {output.directory}
        vamb bin taxvamb {vamb_extra_arg} --taxonomy {input.mmseqs2_out}.formatted --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

