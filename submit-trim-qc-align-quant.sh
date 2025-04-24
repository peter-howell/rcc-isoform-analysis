#!/bin/bash
#SBATCH --job-name=align-and-quantify
#SBATCH --output=logs/nf_%j.out
#SBATCH --error=logs/nf_%j.err
#SBATCH --partition=standard
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1

module load nextflow

export PATH="/home/pjhowell/bin:$PATH"

nextflow -c nextflow.config run trim-align-quant.nf  -resume


