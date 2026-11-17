#!/bin/bash


# Create the combined FASTA file
output_file="working_directory/combined_sequences.fasta"

# Start with an empty output file
> "$output_file"

# Loop through each outlet directory
for sequences_by_outlet in outlet_*; do

    # Get the outlet name
    outlet_name=$(basename "$sequences_by_outlet")

    # Find all FASTA files within that outlet
    find "$sequences_by_outlet" -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) -print0 |
    while IFS= read -r -d '' fasta_file; do

        # Add the outlet name to each FASTA header
        awk -v outlet="$outlet_name" '
        /^>/ {
            sub(/^>/, ">" outlet "_")
        }
        {print}
        ' "$fasta_file" >> "$output_file"

    done
done
