#!/usr/bin/env nextflow

process StringTie {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.gtf")

    script:
    def gtf = file(params.gtf)
    def outdir = file(params.outdir)
    """
    stringtie ${bam} \
        -G ${gtf} \
        -o ${sample_id}.gtf \
        -p ${params.stringtie_cpus} \
        -e
    """
}

Channel
    .fromPath(params.bams)
    .map { bam ->
        def sample_id = bam.baseName
        tuple(sample_id, bam)
    }
    .set { sample_bams }

// Run pipeline steps
workflow {
    sample_bams | StringTie
}

