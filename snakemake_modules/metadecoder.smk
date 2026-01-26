# Metadecoder only takes in sam files: convert the bam files to samfiles
rule bam_to_sam:
    input:
        bamfile = OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam"
    output:
        samfile = temp(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sam.sort"),
    threads: 2
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}" + rulename
    conda: THIS_FILE_DIR / "envs/samtools.yaml"
    shell: 
        """
        samtools view -h {input.bamfile} > {output.samfile} 2> {log}
        """

# Tool only takes decompressed files in
# Additionally the contigs need to be named differently because they cache not to a unique file but to one based on the contig name
# So they will crash internally if not without any error message.
# https://github.com/liu-congcong/MetaDecoder/issues/4
rulename = "decompress_fastafile"
rule decompress_fastafile:
    input:
        contigs = contigs_all,
    output:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
    shell:
        """
        # # Decompress fastafile 
        # gzip --decompress -c {input.contigs} > {output.contigs_decompressed}
        # echo decompression done

        if [[ "{input.contigs}" == *.gz ]]; then
            echo "Decompressing {input.contigs}..."
            gzip --decompress -c {input.contigs} > {output.contigs_decompressed}
        else
            echo "{input.contigs} is not compressed. Copying..."
            cp {input.contigs} {output.contigs_decompressed}
        fi
        """

rulename = "metadecoder_pre"
rule metadecoder_pre:
    input:
        samfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        coverage_file = OUTDIR /  "{key}/metadecoder/coverage_file.coverage",
        seed = OUTDIR /  "{key}/metadecoder/seed.seed",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/meta_decoder.yaml"
    shell:
        """
        # Calculate coverage
        metadecoder coverage --threads {threads} -s {input.samfiles} -o {output.coverage_file}
        echo metadecoder coverage done

        # Get seed
        metadecoder seed --threads {threads} -f {input.contigs_decompressed} -o {output.seed}
        echo metadecoder seed done
        """

rulename = "metadecoder"
rule metadecoder:
    input:
        samfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        coverage_file = OUTDIR /  "{key}/metadecoder/coverage_file.coverage",
        seed = OUTDIR /  "{key}/metadecoder/seed.seed",
    output:
        metadecoder = directory(OUTDIR /  "{key}/metadecoder/clusters/clusters.metadecoder"),
        metadecoder_clusterfile = OUTDIR /  "{key}/metadecoder/clusters/clusters.metadecoder.1.fasta",
        metadecoder_bin_dir = directory(OUTDIR /  "{key}/metadecoder/clusters"),
    threads: 64
    resources: walltime = walltime_fn(rulename), mem_gb = "1000", gpu=""
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/meta_decoder_2.yaml"
    shell:
        """
        # Actually cluster 
        metadecoder cluster -f {input.contigs_decompressed} -c {input.coverage_file} -s {input.seed} -o {output.metadecoder}
        echo metadecoder cluster done
        mkdir -p {output.metadecoder}
        """

rulename = "metadecoder_checkm"
rule metadecoder_checkm:
    input: 
        bin_dir = directory(OUTDIR /  "{key}/metadecoder/clusters"),
        coverage_file = OUTDIR /  "{key}/metadecoder/coverage_file.coverage", # to make sure the rule above is done
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/metadecoder"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fasta' --database_path {params.database}
        """
