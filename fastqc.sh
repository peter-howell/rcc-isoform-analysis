#!/bin/bash
#SBATCH -J FASTQC
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -o logs/fastqc-%j.out
#SBATCH -e logs/fastqc-%j.err

module load fastqc

cd /home/pjhowell/rcc-isoform-analysis
echo "currently in this place:"
pwd

DATADIR="${1:-"data/fastq"}"

FILES=(${DATADIR}/*.fastq)

fastqc -t $SLURM_CPUS_PER_TASK -o qc "${FILES[@]}"


