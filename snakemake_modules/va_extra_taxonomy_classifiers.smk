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
        metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key}.metabuli --tie-ratio 1 --seq-mode 1 --threads {threads}
        # metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key}.metabuli --seq-mode 1 --threads {threads}
        # metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key}.metabuli --accession-level 1 --seq-mode 1 --threads {threads}
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
   # conda: "old_taxconv"
   shell:
       """
       taxconverter metabuli -c {input.metabuli_classification} -r {input.metabuli_report} -o {output.metabuli_classification}
       # python {params.script} {input.metabuli_classification} {input.metabuli_report} > {output.metabuli_classification}
       """



all_targets = {
    "Airways":"/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/new_tax_metabuli/airways_taxonomy_metabuli_otu.tsv",
    "Oral":"/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/new_tax_metabuli/oral_taxonomy_metabuli_otu.tsv",
    "Skin":"/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/new_tax_metabuli/skin_taxonomy_metabuli_otu.tsv",
    "Urogenital":"/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/new_tax_metabuli/urog_taxonomy_metabuli_otu.tsv",
    "Gastrointestinal":"/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/new_tax_metabuli/gi_taxonomy_metabuli_otu.tsv",
}
#
# # Run taxvamb 
# rulename = "run_taxvamb_metabuli"
# rule run_taxvamb_metabuli:
#     input:
#         contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
#         contigs = contigs_all,
#         bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
#         taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
#         # taxonomy =  lambda wildcards: all_targets[wildcards.key]#"{gt_tax}",
#     output:
#         directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
#         out_tax = directory(os.path.join(OUTDIR,"{key}", "metabuli_taxvamb_default_taxonomy_formatted.tsv")),
#         out_id = directory(os.path.join(OUTDIR,"{key}", "id.txt")),
#         out_id_tax = directory(os.path.join(OUTDIR,"{key}", "tax_id.txt")),
#         bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_default','vaevae_clusters_split.tsv'),
#         compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_default/composition.npz'),
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/vamb.yaml"
#     shell:
#         """
#         rm -rf {output.directory}
#
#         # Get list of contig names in each file
#         cat {input.contigs_decompressed} | grep '>' | sed 's/>//' > {output.out_id}
#         cut -f1 {input.taxonomy} > {output.out_id_tax}
#
#         # Format tax
#         echo -e "contigs\tpredictions" > {output.out_tax}
#         cat {input.taxonomy} | awk -F "\t" '{{print $1 "\t" $9 }}' >> {output.out_tax}
#
#         # For missing contigs without tax annotation add tax annotation
#         python /maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/taxvamb_paper_benchmarks/scripts/fix_tax.py {output.out_id} {output.out_id_tax} >> {output.out_tax}
#
#         vamb bin taxvamb --cuda --taxonomy {output.out_tax} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
#         """



# Run taxvamb 
rulename = "run_taxvamb_metabuli"
rule run_taxvamb_metabuli:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_default','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """
# cat metadecoder/Airways_contigs.flt.fna | grep '>' | sed 's/>//' > id.txt


# # Run taxvamb 
#
# rulename = "run_taxvamb_metabuli"
# rule 1_run_taxvamb_metabuli:
#     input:
#         metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_classifications.tsv",
#         metabuli_report= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_report.tsv",
#         contigs = contigs_all,
#         bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
#         # taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
#     output:
#         directory = directory(os.path.join(OUTDIR,"{key}", '1_metabuli_taxvamb_default')),
#         compo = os.path.join(OUTDIR, '{key}','1_metabuli_taxvamb_default/composition.npz'),
#         bins = os.path.join(OUTDIR,"{key}",'1_metabuli_taxvamb_default','vaevae_clusters_split.tsv'),
#         tax = os.path.join(OUTDIR,"{key}",'1_metabuli_taxvamb_default_tax'),
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     # conda: THIS_FILE_DIR / "envs/vamb.yaml"
#     shell:
#         """
#         util format-tax {input.metabuli_report} {input.metabuli_classification} > {output.tax}
#         rm -rf {output.directory}
#         vamb bin taxvamb --cuda --taxonomy {output.tax} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
#         """

# Run taxvamb_no_predictor 
rulename = "run_taxvamb_no_predictor_metabuli"
rule run_taxvamb_no_predictor_metabuli:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """


# rulename = "metabuli_checkm"
# rule metabuli_checkm:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/metabuli_taxvamb"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """
#
# rulename = "metabuli_checkm_no_predictor"
# rule metabuli_checkm_no_predictor:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_no_predictor')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/metabuli_taxvamb_no_predictor"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """
#
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
    # conda: THIS_FILE_DIR / "envs/kraken2.yaml"  
    conda: "/maps/projects/rasmussen/data/taxvamb_benchmarks/taxvamb_benchmarks/airways/.snakemake/conda/28749b7330b451f59560ddf3dac9c4c9_"
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
    # conda: "old_taxconv"
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
        bins = os.path.join(OUTDIR,"{key}",'kraken_taxvamb_default','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','kraken_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        # sed 's/\t$/\tunknown/' {input.taxonomy}  > {input.taxonomy}.fmt
        vamb bin taxvamb --cuda --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """


# Run taxvamb_no_predictor 
rulename = "run_taxvamb_no_predictor_kraken"
rule run_taxvamb_no_predictor_kraken:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/kraken2/taxvamb_formatted_classifications",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'kraken_taxvamb_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','kraken_taxvamb_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# rulename = "kraken_checkm_no_predictor"
# rule kraken_checkm_no_predictor:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_no_predictor')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/kraken_taxvamb_no_predictor"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """
#
# rulename = "kraken_checkm"
# rule kraken_checkm:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_default')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/kraken_taxvamb"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """
#
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
    # conda: "old_taxconv"
    shell:
        """
        taxconverter centrifuge -i {input.centrifuge} -o {output.centrifuge}
        """

# Run taxvamb_no_predictor 
rulename = "run_taxvamb_no_predictor_centrifuge"
rule run_taxvamb_no_predictor_centrifuge:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/centrifuge/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'centrifuge_taxvamb_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','centrifuge_taxvamb_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb_no_predictor 
rulename = "run_taxvamb_centrifuge"
rule run_taxvamb_centrifuge:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/centrifuge/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb')),
        bins = os.path.join(OUTDIR,"{key}",'centrifuge_taxvamb','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','centrifuge_taxvamb/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        # sed 's/\t$/\tunknown/' {input.taxonomy}  > {input.taxonomy}.fmt
        vamb bin taxvamb --cuda  --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb_no_predictor 
rulename = "run_taxvamb_premade"
rule run_taxvamb_premade:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = "/maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/fmt_mmseqs_pred.tsv"
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'premade_taxvamb')),
        bins = os.path.join(OUTDIR,"{key}",'premade_taxvamb','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','premade_taxvamb/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda  --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb 
rulename = "run_taxvamb_gtdb_no_predictor"
rule run_taxvamb_gtdb_no_predictor:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        mmseqs2_out = OUTDIR /  "{key}/classifiers/mmseqs2/gtdb_lca.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'gtdb_taxvamb_default_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'gtdb_taxvamb_default_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','gtdb_taxvamb_default_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        # Format taxonomy file
        echo -e "contigs\tpredictions" > {input.mmseqs2_out}.formatted_nopred
        #sed 's/\t$/\tunknown/' {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted
        cat {input.mmseqs2_out} |  awk -F "\t" '{{print $1 "\t" $9 }}' >> {input.mmseqs2_out}.formatted_nopred
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.mmseqs2_out}.formatted_nopred --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

## Collecting
all_bin_dirs_clas = {
    # "kraken_taxvamb_default": OUTDIR / "{key}/kraken_taxvamb_default/vaevae_clusters_split.tsv",
    # "kraken_taxvamb_no_predictor": OUTDIR / "{key}/kraken_taxvamb_no_predictor/vaevae_clusters_split.tsv",
    "run_taxvamb_centrifuge": OUTDIR / "{key}/centrifuge_taxvamb/vaevae_clusters_split.tsv",
    # "centrifuge_taxvamb_no_predictor": OUTDIR / "{key}/centrifuge_taxvamb_no_predictor/vaevae_clusters_split.tsv",
    # "metabuli_taxvamb_default": OUTDIR / "{key}/metabuli_taxvamb_default/vaevae_clusters_split.tsv",
    # "metabuli_taxvamb_no_predictor": OUTDIR / "{key}/metabuli_taxvamb_no_predictor/vaevae_clusters_split.tsv",
    # "run_taxvamb_gtdb": OUTDIR / "{key}/gtdb_taxvamb_default/vaevae_clusters_split.tsv",
    "run_taxvamb_gtdb_w_unknown": OUTDIR / "{key}/gtdb_taxvamb_default_w_unknown/vaevae_clusters_split.tsv",
    # "run_taxvamb_gtdb_no_predictor": OUTDIR / "{key}/gtdb_taxvamb_default_no_predictor/vaevae_clusters_split.tsv",
    # "kalmari_taxvamb_default": OUTDIR / "{key}/kalmari_taxvamb_default/vaevae_clusters_split.tsv",
    # "trembl_taxvamb_default": OUTDIR / "{key}/trembl_taxvamb_default/vaevae_clusters_split.tsv",
}
bin_dir_names_clas = all_bin_dirs_clas.keys()

# Run taxvamb 
rulename = "format_bins_class"
rule format_bins_class:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        bins = lambda wildcards: all_bin_dirs_clas[wildcards.bins_clas]
    output:
        directory = directory(OUTDIR / "{key}/formatted_vamb_bins/{bins_clas}"),
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
        # expand(OUTDIR / "{{key}}/gunc/{bin_dir}",  bin_dir = bin_dir_names)
    output: 
        OUTDIR / "{key}/tmp/checkm.done"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    shell:
        """
        touch  {output}
        """


rulename = "format_bins_class_filtered"
rule format_bins_class_filtered:
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

rulename = "gunc_2"
rule gunc_2:
    input:
        # semi_bin_decompressed = OUTDIR / "{key}/tmp/semibin_decompressed.done", # recomment me
        bin_dir = directory(OUTDIR / "{key}/formatted_vamb_bins_filtered/{bins_clas}"),
    output: 
        dir = directory(OUTDIR / "{key}/gunc/{bins_clas}"),
        file = OUTDIR / "{key}/gunc/{bins_clas}/GUNC.progenomes_2.1.maxCSS_level.tsv"
    params:
        database = config.get("gunc_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/gunc.yaml"
    shell:
        """
        # Input dir should contain files in fasta format -- need to gzip for some
        rm -rf {output.dir}
        mkdir -p {output.dir}
        gunc run --input_dir {input.bin_dir} --file_suffix .fna --db_file {params.database} --threads {threads} --out_dir {output.dir} 
        """

rulename = "collect_gunc2"
rule collect_gunc2:
    input:
        # expand(OUTDIR / "{{key}}/gunc/{bin_dir}",  bin_dir = bin_dir_names)
        expand(OUTDIR / "{{key}}/gunc/{bins_clas}/GUNC.progenomes_2.1.maxCSS_level.tsv", bins_clas = bin_dir_names_clas)
    output: 
        OUTDIR / "{key}/tmp/gunc_2.done"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    shell:
        """
        touch  {output}
        """

# rulename = "centrifuge_checkm"
# rule centrifuge_checkm:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_default')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """


# rulename = "centrifuge_checkm_no_predictor"
# rule centrifuge_checkm_no_predictor:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_no_predictor')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb_no_predictor"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """
#
# rulename = "centrifuge_checkm"
# rule centrifuge_checkm:
#     input: 
#         bin_dir = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_default')),
#     output:
#         outdir = directory(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb"),
#     threads: threads_fn(rulename)
#     params:
#         database = config.get("checkm2_database")
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/checkm2.yaml"
#     shell:
#         """
#         checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
#         """



## Collecting
# all_bin_dirs_clas_2 = {
#     # "metabuli_taxvamb_default": OUTDIR / "{key}/metabuli_taxvamb_default/vaevae_clusters_split.tsv",
#     # "metabuli_taxvamb_no_predictor": OUTDIR / "{key}/metabuli_taxvamb_no_predictor/vaevae_clusters_split.tsv",
#     # "run_taxvamb_premade": OUTDIR / "{key}/premade_taxvamb/vaevae_clusters_split.tsv",
#     # "run_taxvamb_gtdb": OUTDIR / "{key}/gtdb_taxvamb_default/vaevae_clusters_split.tsv",
#     # "run_taxvamb_gtdb_no_predictor": OUTDIR / "{key}/gtdb_taxvamb_default_no_predictor/vaevae_clusters_split.tsv",
#     # "run_taxvamb_centrifuge": OUTDIR / "{key}/centrifuge_taxvamb/vaevae_clusters_split.tsv",
#     # "kraken_taxvamb_default": OUTDIR / "{key}/kraken_taxvamb_default/vaevae_clusters_split.tsv",
#     "kraken_taxvamb_no_predictor": OUTDIR / "{key}/kraken_taxvamb_no_predictor/vaevae_clusters_split.tsv",
#     "centrifuge_taxvamb_no_predictor": OUTDIR / "{key}/centrifuge_taxvamb_no_predictor/vaevae_clusters_split.tsv",
# }
#
# rulename = "run_all"
# rule run_all:
#     input:
#         expand("{bins_clas_2}", bins_clas_2 = all_bin_dirs_clas_2.values())
#         # expand(OUTDIR / "{{key}}/gunc/{bin_dir}",  bin_dir = bin_dir_names)
#     output: 
#         OUTDIR / "{key}/tmp/run_all.done"
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     shell:
#         """
#         touch  {output}
#         """
