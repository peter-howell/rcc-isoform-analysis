#!/bin/bash
#SBATCH --job-name=align-and-quantify
#SBATCH --output=logs/nf_%j.out
#SBATCH --error=logs/nf_%j.err
#SBATCH --partition=standard
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load nextflow
nextflow run main.nf -profile slurm -resume


