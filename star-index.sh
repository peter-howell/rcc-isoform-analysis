#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -n 1
#SBATCH -o logs/star-index-%j.out
#SBATCH -e logs/star-index-%j.err
#SBATCH -J STAR-INDEX

export PATH="/home/pjhowell/bin:$PATH"

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles data/ref/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile data/ref/gencode.v47.primary_assembly.annotation.gtf \
     --sjdbOverhang 100

