#!/usr/bin/env bash

minimap2 \
    -t {threads} \
    --secondary=no \
    -I 8G --eqx --MD -Y -y \
    -ax lr:hq \
    {input.dsa} \
    <(samtools fastq -@ 32 -T "*" $INBAM) |
    samtools view -@ 32 -u \
        >{output.bam}
