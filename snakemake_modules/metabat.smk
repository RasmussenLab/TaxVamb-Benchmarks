# Contigs are given outside of the scope of the output directory
# metabat uses singularity which mounts only the output directory
# we therefore need to copy the contigs inside the output directory to be mounted properbly
rulename = "tmp_copy"
rule tmp_copy:
    input:
        contigs = contigs_all,
    output: 
        tmp_contigs = temp("{key}/contigs.fasta.gz"),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    shell:
        """
        cp {input.contigs} {output.tmp_contigs}
        """

rulename = "metabat"
rule metabat:
    input:
        contigs = "{key}/contigs.fasta.gz",
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
    output: 
        depht = OUTDIR /  "{key}/metabat/depht.txt",
        metabat = directory(OUTDIR /  "{key}/metabat/metabat"),
    threads: threads_fn(rulename)
    params:
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    container: "docker://metabat/metabat:v2.17-66-ga512006"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depht} {input.bamfiles}
        mkdir -p {output.metabat}
        metabat2 -i {input.contigs} -a {output.depht} -o {output.metabat}/metabat
        """


rulename = "metabat_checkm"
rule metabat_checkm:
    input:
        bin_dir = OUTDIR /  "{key}/metabat/metabat",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/metabat"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: BENCHMARK_DIR + "{key}_" + rulename
    log: LOG_DIR + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension '.fa' --database_path {params.database}
        """





#     (base) [bxc755@esrumcmpn05fl taxvamb_benchmarks]$ /maps/projects/rasmussen/people/bxc755/conda_env/conda/bin/snakemake --snakefile /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/taxvamb_paper_benchmarks/snakefile.smk --rerun-triggers mtime --nolock -c 40 -p --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --keep-incomplete --config bam_contig=../taxvamb_paper_benchmarks/data_configs/esrum_airways_bamfile_contigs.tsv output_directory=airways --direct ory airways

    # metabat2 --help
    # runMetaBat.sh --numThreads {threads} --inFile {input.contigs} --outFile {output.metabat} {input.bamfiles} 
# MetaBAT: Metagenome Binning based on Abundance and Tetranucleotide frequency (version 2:2.17.66-ga512006-dirty; 20250124_080221)
# by Don Kang (ddkang@lbl.gov), Feng Li, Jeff Froula, Rob Egan, and Zhong Wang (zhongwang@lbl.gov)
#
# Allowed options:
#   -h [ --help ]                     produce help message
#   -i [ --inFile ] arg               Contigs in (gzipped) fasta file format [Mandatory]
#   -o [ --outFile ] arg              Base file name and path for each bin. The default output is fasta format.
#                                     Use -l option to output only contig names [Mandatory].
#   -a [ --abdFile ] arg              A file having mean and variance of base coverage depth (tab delimited;
#                                     the first column should be contig names, and the first row will be
#                                     considered as the header and be skipped) [Optional].
#   -m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
#   --minSmallContig arg (=1000)      Minimum size of a small contig for recruiting into established bins
#                                     (should be >=500)
#   --maxP arg (=95)                  Percentage of 'good' contigs considered for binning decided by connection
#                                     among contigs. The greater, the more sensitive.
#   --minS arg (=60)                  Minimum score of a edge for binning (should be between 1 and 99). The
#                                     greater, the more specific.
#   --maxEdges arg (=200)             Maximum number of edges per node. The greater, the more sensitive.
#   --pTNF arg (=0)                   TNF probability cutoff for building TNF graph. Use a %% value between 1
#                                     and 100 to skip the auto preparation step. (0: auto).
#   --noAdd                           Turning off additional binning for lost or small contigs.
#   --minRecruitingSize arg (=10)     (if not noAdd) Minimum cluster size for recruiting of small and leftover
#                                     contigs
#   --recruitToAbdCentroid            [EXPERIMENTAL] If set (and not noAdd), use the weighted-by-abundance
#                                     centroid of a cluster to recruit small and lost contigs.  Potentially
#                                     reduces sensitivity and improves speed
#   --recruitWithTNF arg (=0)         [EXPERIMENTAL] If non-zero (and not noAdd), uses this factor against the
#                                     large-contig TNF threshold to require small and lost contigs have at
#                                     least that TNF distance from the centroid of the recruiting cluster.
#                                     Recommend 0.9-1.0.
#   --cvExt                           When a coverage file without variance (from third party tools) is used
#                                     instead of abdFile from jgi_summarize_bam_contig_depths.
#   -x [ --minCV ] arg (=1)           Minimum mean coverage of a contig in each library for binning.
#   --minCVSum arg (=1)               Minimum total effective mean coverage of a contig (sum of depth over
#                                     minCV) for binning.
#   -s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
#   -t [ --numThreads ] arg (=0)      Number of threads to use (0: use all cores).
#   -l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences.
#   --saveCls                         Save cluster memberships as a matrix format
#   --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
#   --noBinOut                        No bin output. Usually combined with --saveCls to check only contig
#                                     memberships
#   --noSampleDepths                  Do not include per-sample depths in bin fasta headers
#   --seed arg (=0)                   For exact reproducibility. (0: use random seed)
#   -d [ --debug ]                    Debug output
#   -q [ --quiet ]                    Be less verbose verbose output
#   -v [ --verbose ]                  Be more verbose in output (on by default)
