#!/bin/bash

# Declare an associative array to store the newest file per name
declare -A newest_files

# Find all matching files
while IFS= read -r -d '' file; do
    filename=$(basename "$file")

    # If we haven't seen this filename before, or it's newer than the stored one
    if [[ ! -v newest_files["$filename"] ]] || [[ "$file" -nt "${newest_files[$filename]}" ]]; then
        newest_files["$filename"]="$file"
    fi
done < <(find . -type f -name 'SRR*.Aligned.sortedByCoord.out.bam' -print0)

# Output the selected files
echo "Newest matching files:"
for file in "${newest_files[@]}"; do
    echo "$file"
done

