#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq = 'data/fastq/*.fastq.gz'
params.outdir = 'results'
params.star_index = 'star_index'
params.gtf = 'annotation.gtf'
params.genome = 'genome.fa'
params.threads = 8

process FastQC {
    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.html"), path("*.zip") into fastqc_results

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads.join(' ')}
    """
}

process TrimGalore {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz") into trimmed_reads

    script:
    """
    trim_galore --paired ${reads[0]} ${reads[1]} --gzip -o .
    """
}

process STARAlign {
    input:
    tuple val(sample_id), path(reads)
    path params.star_index
    path params.gtf

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam") into aligned_bams

    script:
    """
    STAR \
        --genomeDir ${params.star_index} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --runThreadN ${params.threads} \
        --outFileNamePrefix ${sample_id}.
    """
}

process StringTie {
    input:
    tuple val(sample_id), path(bam)
    path params.gtf

    output:
    tuple val(sample_id), path("${sample_id}.gtf") into stringtie_gtfs

    script:
    """
    stringtie ${bam} \
        -G ${params.gtf} \
        -o ${sample_id}.gtf \
        -p ${params.threads}
    """
}

// Channel to group reads into pairs
Channel
    .fromFilePairs(params.fastq, size: 2)
    .set { read_pairs }

// Run pipeline steps
workflow {
    read_pairs |>
        FastQC |>
        TrimGalore |>
        STARAlign |>
        StringTie
}

