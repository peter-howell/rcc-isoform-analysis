#!/bin/bash


filename="${1:-"SRR_Acc_List.txt"}"

echo "using accessions listed in $filename"

module load sratoolkit


module load parallel


