#!/bin/bash

#Loop through each outlet directory
for sequences_by_outlet in outlet_*; do
	
	#get the outlet name
	outlet_name=$(basename "$sequences_by_outlet")

	# Find all FASTA files within that outlet
    	find "$sequences_by_outlet" -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) -print0 |
    	while IFS= read -r -d '' fasta_file; do

        # Create a unique output filename
        relative_path="${fasta_file#"$sequences_by_outlet"/}"
        safe_name="${relative_path//\//_}"
        output_file="working_directory/${outlet_name}_${safe_name}"

        # Add the outlet name to each FASTA header
        awk -v outlet="$outlet_name" '
        /^>/ {
            sub(/^>/, ">" outlet "_")
        }
        {print}
        ' "$fasta_file" > "$output_file"

    done
done
