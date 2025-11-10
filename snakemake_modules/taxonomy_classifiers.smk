
# rulename = "metabuli"
# rule metabuli:
#     input:
#         contigs = contigs_all,
#     output:
#         metabuli= OUTDIR /  "{key}/classifiers/metabuli",
#     params: 
#         database = config.get("database")
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/metabuli.yaml"
#     shell:
#         """
#         metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key} --seq-mode 1 --threads {threads}
#         """
# NOTE: does it work with gzipped fastafiles? 

# metabuli classify <i:FASTA/Q> <i:DBDIR> <o:OUTDIR> <Job ID> [options]
# - INPUT : FASTA/Q file of reads you want to classify. (gzip supported)
# - DBDIR : The directory of reference DB. 
# - OUTDIR : The directory where the result files will be generated.
# - Job ID: It will be the prefix of result files.  
#
# # Paired-end
# metabuli classify read_1.fna read_2.fna dbdir outdir jobid
#
# # Single-end
# metabuli classify --seq-mode 1 read.fna dbdir outdir jobid
#
# # Long-read 
# metabuli classify --seq-mode 3 read.fna dbdir outdir jobid
#
#   * Important parameters:
#    --threads : The number of threads used (all by default)
#    --max-ram : The maximum RAM usage. (128 GiB by default)
#    --min-score : The minimum score to be classified 
#    --min-sp-score : The minimum score to be classified at or below species rank. 
#    --taxonomy-path: Directory where the taxonomy dump files are stored. (DBDIR/taxonomy by default)
#    --accession-level : Set 1 to use accession level classification (0 by default). 
#                        It is available when the DB is also built with accession level taxonomy.

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
        bins = os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default','vaevae_clusters_split.tsv'),
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
        # Format taxonomy file
        # echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        # awk -F "\t" '{{print $1 "\t" $5 }}' {input.mmseqs2_out} >> {input.mmseqs2_out}.formatted
        # sed 's/\t$/\tunknown/' {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $5 }}' >> {input.mmseqs2_out}.formatted

        cut -f1,5 {input.mmseqs2_out} > {input.mmseqs2_out}.cut.tsv
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        python {params.script} {input.mmseqs2_out}.cut.tsv >> {input.mmseqs2_out}.formatted

        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.mmseqs2_out}.formatted --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
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
        bins = os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default_w_unknown','vaevae_clusters_split.tsv'),
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
        sed 's/\t$/\tunknown/' {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted_w_unknown
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.mmseqs2_out}.formatted_w_unknown --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb 
rulename = "run_taxvamb_gtdb"
rule run_taxvamb_gtdb:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/gtdb_lca.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','gtdb_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        # Format taxonomy file
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        # sed 's/\t$/\tunknown/' {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted
        cat {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.mmseqs2_out}.formatted --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
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
        bins = os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','trembl_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    params: script =THIS_FILE_DIR / "files_used_in_snakemake_workflow/format_trembl_kalmari.py"
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        # Format taxonomy file
        # echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        # sed 's/\t$/\tunknown/' {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted
        # awk -F "\t" '{{print $1 "\t" $9 }}' {input.mmseqs2_out} >> {input.mmseqs2_out}.formatted

        cut -f1,9 {input.mmseqs2_out} > {input.mmseqs2_out}.cut.tsv
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted
        python {params.script} {input.mmseqs2_out}.cut.tsv >> {input.mmseqs2_out}.formatted

        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.mmseqs2_out}.formatted --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb 
rulename = "format_bins_taxvamb_kalmari"
rule format_bins_taxvamb_kalmari:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'kalmari_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default','vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kalmari_taxvamb_default_bins')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """


# Run taxvamb 
rulename = "format_bins_taxvamb_filtered_kalmari"
rule format_bins_taxvamb_filtered_kalmari:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'kalmari_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'kalmari_taxvamb_default','vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kalmari_taxvamb_default_bins_filtered')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """

# Run taxvamb 
rulename = "format_bins_taxvamb_gtdb"
rule format_bins_taxvamb_gtdb:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default','vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default_bins')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """


# Run taxvamb 
rulename = "format_bins_taxvamb_filtered_gtdb"
rule format_bins_taxvamb_filtered_gtdb:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default','vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default_bins_filtered')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """

# Run taxvamb 
rulename = "format_bins_taxvamb_trembl"
rule format_bins_taxvamb_trembl:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'trembl_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default','vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'trembl_taxvamb_default_bins')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 0 {output.directory} 
        """


# Run taxvamb 
rulename = "format_bins_taxvamb_filtered_trembl"
rule format_bins_taxvamb_filtered_trembl:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'trembl_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'trembl_taxvamb_default','vaevae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'trembl_taxvamb_default_bins_filtered')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 200000 {output.directory} 
        """


rulename = "trembl_checkm"
rule trembl_checkm:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'trembl_taxvamb_default_bins')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/trembl_taxvamb"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """


rulename = "kalmari_checkm"
rule kalmari_checkm:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'kalmari_taxvamb_default_bins')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/kalmari_taxvamb"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """


rulename = "gtdb_checkm"
rule gtdb_checkm:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default_bins')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/gtdb_taxvamb"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

# awk -F "\t" '{print $1 "\t" $5 }' kalmari_lca.tsv


# mmseqs_db_gtdb: "/
# mmseqs_db_trembl: #
# mmseqs_db_kalmari:# rulename = "centrifuge"
# rule centrifuge:
#     output:
#         centrifuge = OUTDIR /  "{key}/classifiers/centrifuge",
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/centrifuge.yaml"
#     shell:
#         """
#
#         -k 1
#         """
#
# rulename = "kraken2"
# rule kraken2:
#     output:
#         kraken2 = OUTDIR /  "{key}/classifiers/kraken2",
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/kraken2.yaml"
#     shell:
#         """
#         echo
#         """
#
