rule pbmm2:
    input:
        dsa=DSA,
        bam=BAM,
    output:
        bam=temp("temp/{sm}.bam"),
        bai=temp("temp/{sm}.bam.bai"),
    conda:
        DEFAULT_ENV
    resources:
        runtime=12 * 60,
        mem_mb=3 * MAX_THREADS * 1024,
        tmpdir="temp/tmp.pbmm2.sorts/",
    threads: MAX_THREADS
    params:
        sort_threads=MAX_THREADS // 4,
        sort_memory=8,
        chunk_size=MAX_THREADS * 3,
    shell:
        "mkdir -p {resources.tmpdir} && "
        " pbmm2 align -j {threads} -J {params.sort_threads}"
        " --preset CCS --sort"
        " --sort-memory {params.sort_memory}G"
        " --chunk-size {params.chunk-size}"
        " --log-level DEBUG"
        " --strip"
        " --unmapped"
        " {input.dsa} {input.bam} {output.bam}"


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
        min_mapq=config.get("min_mapq", 1),
        script=workflow.source_path("../scripts/haplotag-reads-by-asm.py"),
    shell:
        "python {params.script} {input.bam} - {output.assignments}"
        " -t {threads} --hap1-tag {H1_TAG} --hap2-tag {H2_TAG} -m {params.min_mapq}"
        " | samtools view -C -@ {threads} -T {input.dsa}"
        "  --output-fmt-option embed_ref=1"
        "  --write-index -o {output.cram}"


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
        "python {params.script} -t {threads} {input.bam} {input.assignments} -"
        " | samtools view -C -@ {threads} -T {input.ref}"
        "  --output-fmt-option embed_ref=1"
        "  --write-index -o {output.cram}"


rule haplotagged_vcf:
    input:
        vcf=VCF,
        ref=REF,
        cram=rules.add_assignments.output.cram,
        crai=rules.add_assignments.output.crai,
    output:
        vcf="results/{sm}.reference.vcf.gz",
    conda:
        DEFAULT_ENV
    threads: 4
    resources:
        mem_mb=32 * 1024,
        runtime=8 * 60,
    benchmark:
        "benchmark/{sm}/haplotagged_vcf/{sm}.bench.txt"
    shell:
        "whatshap haplotagphase -r {input.ref}"
        " -o {output.vcf} {input.vcf} {input.cram}"
