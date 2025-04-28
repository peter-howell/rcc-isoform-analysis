#!/usr/bin/env nextflow

process TrimGalore {
    tag "$sample_id"

    // cpus params.threads

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz") 

    script:
    """
    trim_galore --cores ${params.trim_cpus} --paired ${reads[0]} ${reads[1]} --gzip -o .
    """
}

process STARAlign {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam") 

    script:
    def index_dir = file(params.star_index)
    """
    STAR \
        --genomeDir ${index_dir} \
        --readFilesIn ${reads1} ${reads2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --runThreadN ${params.star_cpus} \
        --outFileNamePrefix ${sample_id}.
    """
}

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
        -p ${params.stringtie_cpus}
    """
}

// Channel to group reads into pairs
Channel
    .fromFilePairs(params.fastq, size: 2)
    .set { read_pairs }


// Run pipeline steps
workflow {
    read_pairs | TrimGalore | STARAlign | StringTie
}

