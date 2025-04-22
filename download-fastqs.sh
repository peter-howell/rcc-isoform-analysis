#!/bin/bash


filename="${1:-"SRR_Acc_List.txt"}"

echo "using accessions listed in $filename"

module load sratoolkit


while IFS= read -r sra_id; do
	echo "Starting on accession $sra_id"
	prefetch -p "$sra_id"
	echo "done prefetching $sra_id. beginning fasterq-dump"
	fasterq-dump -O data/fastq -e 2 "$sra_id" 

	echo "completed $sra_id"
done < "$filename"


echo "all done"

