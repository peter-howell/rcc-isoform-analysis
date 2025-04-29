#!/bin/bash




FILE="${1:-"bams.txt"}"




while IFS= read -r file_name; do
	base_file_name="$(basename "$file_name")"
	cp "work/$file_name" "results/bam/$base_file_name"
done < "$FILE"

echo "ALL DONE"
