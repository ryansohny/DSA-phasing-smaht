# Process each input file individually
rule extract_fastq:
    input:
        bam=get_bam,
    output:
        fastq=pipe("temp/{sm}.{file_idx}.fastq"),
    conda:
        DEFAULT_ENV
    threads: 8
    shell:
        'samtools fastq -@ {threads} -T "*" {input.bam} > {output.fastq}'


rule align:
    input:
        dsa=get_dsa,
        fastq=rules.extract_fastq.output.fastq,
        bam=get_bam,
    output:
        bam=pipe("temp/{sm}.{file_idx}.bam"),
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=MAX_THREADS * 1024,
    params:
        sample=bam_header_sm_settings,
        mm2_preset=config.get("mm2_preset", "'lr:hq'"),
        mm2_extra_opts=config.get("mm2_extra_options", ""),
    threads: MAX_THREADS
    shell:
        "minimap2"
        " -t {threads}"
        " --secondary=no -I 8G --eqx --MD -Y -y"
        " -ax {params.mm2_preset}"
        " {params.mm2_extra_opts}"
        " {input.dsa} {input.fastq}"
        " | rb add-rg -u {params.sample} -t {threads} {input.bam}"
        " > {output.bam}"


rule haplotag_and_sort:
    input:
        dsa=get_dsa,
        bam=rules.align.output.bam,
    output:
        assignments=read_assignment_result(),
        cram=temp("temp/{sm}.{file_idx}.dsa.cram"),
        crai=temp("temp/{sm}.{file_idx}.dsa.cram.crai"),
    conda:
        DEFAULT_ENV
    threads: 8
    resources:
        mem_mb=(8 * 8 + 8) * 1024,
    params:
        min_mapq=config.get("min_mapq", 1),
        script=workflow.source_path("../scripts/haplotag-reads-by-asm.py"),
        sort_memory=8,  # GB per thread
        h1_tag=get_h1_tag,
        h2_tag=get_h2_tag,
    shell:
        "python {params.script} {input.bam} - {output.assignments}"
        " -t {threads} -m {params.min_mapq}"
        " --hap1-tag {params.h1_tag} --hap2-tag {params.h2_tag}"
        " | samtools sort -u -@ {threads} -m {params.sort_memory}G"
        " | samtools view -C -@ {threads} -T {input.dsa}"
        "  --output-fmt-option embed_ref=1"
        "  --write-index -o {output.cram}"


# For ONT: run modkit on each file individually
rule modkit:
    input:
        cram=rules.haplotag_and_sort.output.cram,
        dsa=get_dsa,
    output:
        aln=temp("temp/{sm}.{file_idx}.modkit.dsa.bam"),
    log:
        "results/modkit/{sm}.{file_idx}.modkit.log",
    conda:
        DEFAULT_ENV
    resources:
        runtime=24 * 60,
        mem_mb=16 * 1024,
    threads: 8
    params:
        ft_nuc=config.get("ft_nuc_params", ""),
    shell:
        "modkit call-mods --log-filepath {log} -t {threads} -p 0.1 {input.cram} -"
        " | ft add-nucleosomes -t {threads} {params.ft_nuc} - {output.aln}"


rule fire:
    input:
        cram=get_fire_input,
        dsa=get_dsa,
    output:
        cram=temp("temp/{sm}.{file_idx}.fire.dsa.cram"),
        crai=temp("temp/{sm}.{file_idx}.fire.dsa.cram.crai"),
    conda:
        DEFAULT_ENV
    threads: 16
    resources:
        runtime=12 * 60,
        mem_mb=16 * 1024,
    params:
        use_ont=lambda wc: "--ont" if is_ont(wc) else "",
    shell:
        """
        ft fire -u -t {threads} {params.use_ont} {input.cram} \
            | samtools view -C -@ {threads} -T {input.dsa} \
                --output-fmt-option embed_ref=1 \
                --write-index -o {output.cram}
        """


# Merge all files for a sample (after FIRE)
rule merge_sample:
    input:
        crams=get_crams_to_merge,
        crais=get_crais_to_merge,
        dsa=get_dsa,
    output:
        cram="results/{sm}.dsa.cram",
        crai="results/{sm}.dsa.cram.crai",
    conda:
        DEFAULT_ENV
    threads: 16
    resources:
        mem_mb=16 * 1024,
    run:
        # If only one file, just copy it instead of merging
        if len(input.crams) == 1:
            shell("cp {input.crams[0]} {output.cram}")
            shell("cp {input.crais[0]} {output.crai}")
        else:
            # Multiple files: merge them
            shell(
                "samtools merge -@ {threads} --write-index"
                " --reference {input.dsa}"
                " --output-fmt cram"
                " --output-fmt-option embed_ref=1"
                " -o {output.cram}"
                " {input.crams}"
            )


# add read assignments to the input bam files
rule qc:
    input:
        cram=get_final_cram,
    output:
        txt="results/{sm}.qc.tbl.gz",
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=16 * 1024,
    threads: 16
    shell:
        "ft qc --acf -t {threads} {input.cram} {output.txt}"


# realign to shared reference if provided
rule realign_to_shared_ref:
    input:
        ref=SHARED_REF,
        bam=get_final_cram,
    output:
        cram="results/{sm}.shared.ref.cram",
        crai="results/{sm}.shared.ref.cram.crai",
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=(MAX_THREADS * 4 + 8) * 1024,
    params:
        sort_memory=3,  # GB per thread
        sample=bam_header_sm_settings,
        mm2_preset=config.get("mm2_preset", "'lr:hq'"),
        mm2_extra_opts=config.get("mm2_extra_options", ""),
    threads: MAX_THREADS
    shell:
        "minimap2"
        " -t {threads}"
        " --secondary=no -I 8G --eqx --MD -Y -y"
        " -ax {params.mm2_preset}"
        " {params.mm2_extra_opts}"
        ' {input.ref} <(samtools fastq -@ {threads} -T "*" {input.bam})'
        " | rb add-rg -u {params.sample} -t {threads} {input.bam}"
        " | samtools sort -u -@ {threads} -m {params.sort_memory}G"
        " | samtools view -C -@ {threads} -T {input.ref}"
        "  --output-fmt-option embed_ref=1"
        "  --write-index -o {output.cram}"
