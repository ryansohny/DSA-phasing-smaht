rule align:
    input:
        dsa=get_dsa,
        bam=get_bam,
    output:
        bam=pipe("temp/{sm}.bam"),
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=MAX_THREADS * 1024,
    threads: MAX_THREADS
    shell:
        "minimap2"
        " -t {threads}"
        " --secondary=no -I 8G --eqx --MD -Y -y"
        " -ax 'lr:hq'"
        ' {input.dsa} <(samtools fastq -@ {threads} -T "*" {input.bam})'
        " | samtools view -u -@ {threads} > {output.bam}"


rule haplotag_and_sort:
    input:
        dsa=get_dsa,
        bam=rules.align.output.bam,
    output:
        assignments="results/{sm}.assignments.tsv.gz",
        cram="results/{sm}.dsa.cram",
        crai="results/{sm}.dsa.cram.crai",
    conda:
        DEFAULT_ENV
    threads: MAX_THREADS // 4
    params:
        min_mapq=config.get("min_mapq", 1),
        script=workflow.source_path("../scripts/haplotag-reads-by-asm.py"),
        mem_mb=(MAX_THREADS * 4 + 8) * 1024,
        sort_memory=4,  # GB per thread
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

rule modkit:
    input:
        cram=rules.haplotag_and_sort.output.cram,
        dsa=get_dsa,
    output:
        bam="results/{sm}.modkit.dsa.cram",
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=16 * 1024,
    threads: 16
    params:
        ft_nuc=config.get("ft_nuc_params", ""),
    shell:
        "modkit call-mods -t {threads} -p 0.1 {input.cram} -"
        " | ft add-nucleosomes -t {threads} {params.ft_nuc}"
        " | samtools view -C -@ {threads} -T {input.dsa}"
        " --output-fmt-option embed_ref=1"
        " --write-index -o {output.bam}"


# add read assignments to the input bam files
rule phase_percentages:
    input:
        cram=rules.haplotag_and_sort.output.cram,
    output:
        txt="results/{sm}.phase-percentages.txt",
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=16 * 1024,
    threads: 16
    params:
        script=workflow.source_path("../scripts/phase-percentages.py"),
    shell:
        "python {params.script} -t {threads} {input.cram} > {output.txt}"
