rule haplotagged_vcf:
    input:
        bam=rules.haplotagged_bam.output.bam,
        bai=rules.haplotagged_bai.output.bai,
        vcf=rules.hiphase_vcf.output.vcf,
        ref=get_ref,
    output:
        vcf="results/{sm}/{sm}.haplotagged.vcf.gz",
    conda:
        CONDA
    threads: 4
    resources:
        mem_mb=32 * 1024,
        runtime=8 * 60,
    benchmark:
        "benchmark/{sm}/haplotagged_vcf/{sm}.bench.txt"
    shell:
        """
        whatshap haplotagphase \
            -r {input.ref} \
            -o {output.vcf} \
            {input.vcf} {input.bam}
        """
