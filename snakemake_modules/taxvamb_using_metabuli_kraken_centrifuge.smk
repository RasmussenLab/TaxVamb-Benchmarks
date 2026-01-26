########### METABULI ####################

rulename = "metabuli"
rule metabuli:
    input:
        contigs = contigs_all,
    output:
        metabuli= directory(OUTDIR /  "{key}/classifiers/metabuli"),
        metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_classifications.tsv",
        metabuli_report= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_report.tsv",
    params: 
        database = config.get("metabuli_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/metabuli.yaml"
    shell:
        """
        metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key}.metabuli --seq-mode 1 --threads {threads}
        """

rulename = "metabuli_taxconv"
rule metabuli_taxconv:
   input:
       metabuli= directory(OUTDIR /  "{key}/classifiers/metabuli"),
       metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_classifications.tsv",
       metabuli_report= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_report.tsv",
   output:
       metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
   threads: threads_fn(rulename)
   params: script =THIS_FILE_DIR / "files_used_in_snakemake_workflow/format_metabuli.py"
   resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
   benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
   log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
   conda: "taxconv"
   shell:
       """
       taxconverter metabuli -c {input.metabuli_classification} -r {input.metabuli_report} -o {output.metabuli_classification}
       """

# Run taxvamb 
rulename = "run_taxvamb_metabuli"
rule run_taxvamb_metabuli:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
        bins = OUTDIR / "{key}/metabuli_taxvamb_default/vaevae_clusters_split.tsv",
        compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb {vamb_extra_arg} --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

########### KRAKEN2 ####################

rulename = "kraken2"
rule kraken2:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2/kraken2_predictions.tsv",
    params: 
        database = config.get("kraken2_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/kraken2.yaml"  
    shell:
        """
        kraken2 --minimum-hit-groups 3 --db {params.database} --threads {threads}  {input.contigs_decompressed} > {output.kraken2}
        """

rulename = "kraken_taxconv"
rule kraken_taxconv:
    input:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2/kraken2_predictions.tsv",
    output:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2/taxvamb_formatted_classifications",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: "taxconv"
    shell:
        """
        taxconverter kraken2 -i {input.kraken2} -o {output.kraken2}
        """

# Run taxvamb 
rulename = "run_taxvamb_kraken"
rule run_taxvamb_kraken:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/kraken2/taxvamb_formatted_classifications",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_default')),
        bins = OUTDIR / "{key}/kraken_taxvamb_default/vaevae_clusters_split.tsv",
        compo = os.path.join(OUTDIR, '{key}','kraken_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        # sed 's/\t$/\tunknown/' {input.taxonomy}  > {input.taxonomy}.fmt  # Taxonomy levels not classified should output 'unknown' instaed of an empty string ('')
        vamb bin taxvamb {vamb_extra_arg} --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

########### CENTRIFUGE ####################

rulename = "centrifuge"
rule centrifuge:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge/centrifuge_predictions.tsv",
    params: 
        database_dir = config.get("centrifuge_database_dir"),
        database_name = config.get("centrifuge_database_name")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/centrifuge.yaml"
    shell:
        """
        centrifuge -x {params.database_dir}/{params.database_name} -k 1 -f {input.contigs_decompressed} --threads {threads} > {output.centrifuge}
        """

rulename = "centri_taxconv"
rule centri_taxconv:
    input:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge/centrifuge_predictions.tsv",
    output:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge/taxvamb_formatted_classifications.tsv",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: "taxconv"
    shell:
        """
        taxconverter centrifuge -i {input.centrifuge} -o {output.centrifuge}
        """


rulename = "run_taxvamb_centrifuge"
rule run_taxvamb_centrifuge:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/centrifuge/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb')),
        bins = OUTDIR / "{key}/centrifuge_taxvamb/vaevae_clusters_split.tsv",
        compo = os.path.join(OUTDIR, '{key}','centrifuge_taxvamb/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        # sed 's/\t$/\tunknown/' {input.taxonomy}  > {input.taxonomy}.fmt        # Taxonomy levels not classified should output 'unknown' instaed of an empty string ('')
        vamb bin taxvamb {vamb_extra_arg}  --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

