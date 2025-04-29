#!/bin/bash


#SBATCH -J quant
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -o logs/prepDE-%j.out



export PATH="/home/pjhowell/bin:$PATH"

source /home/pjhowell/miniconda3/bin/activate base

python /home/pjhowell/src/stringtie-3.0.0/prepDE.py3 DE-samples.csv
mv transcript_count_matrix.csv results/
mv gene_count_matrix.csv results/
echo "ALL DONE"

exit 0

