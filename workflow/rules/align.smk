rule pbmm2:
    input:
        dsa=DSA,
        bam=BAM,
    output:
        bam=temp("temp/{sm}.bam"),
        bai=temp("temp/{sm}.bam.bai"),
        pbi=temp("temp/{sm}.bam.pbi"),
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=2 * MAX_THREADS * 1024,
    threads: MAX_THREADS
    shell:
        """
        pbmm2 align \
            -j {threads} \
            --preset CCS --sort \
            --sort-memory 1G \
            --log-level INFO \
            --strip \
            --unmapped \
            {input.ref} {input.bam} {output.bam} 
        """

# --sample "{wildcards.sm}" \


# convert to cram
rule cram:
    input:
        dsa=DSA,
        bam=rules.pbmm2.output.bam,
    output:
        assignments="results/{sm}.assignments.tsv.gz",
        cram="results/{sm}.dsa.cram",
        crai="results/{sm}.dsa.cram.crai",
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=16 * 1024,
    threads: 32
    params:
        script=workflow.source_path("../scripts/haplotag-reads-by-asm.py"),
    shell:
        "python {params.script} {input.bam} - {output.assignments}"
        " -t {threads} --hap1-tag {H1_TAG} --hap2-tag {H2_TAG}"
        "| samtools view -C -@ {threads} -T {input.dsa}"
        " --output-fmt-option embed_ref=1"
        " --write-index -o {output.cram}"


# add read assignments to the input bam files
rule add_assignments:
    input:
        assignments=rules.cram.output.assignments,
        bam=BAM,
        ref=REF,
    output:
        cram="results/{sm}.reference.cram",
        crai="results/{sm}.reference.cram.crai",
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=16 * 1024,
    threads: 32
    params:
        script=workflow.source_path("../scripts/add-assignments-to-bam.py"),
    shell:
        "python {params.script} {input.bam} {input.assignments} -"
        "| samtools view -C -@ {threads} -T {input.ref}"
        " --output-fmt-option embed_ref=1"
        " --write-index -o {output.cram}"

