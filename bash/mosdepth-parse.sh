#!/bin/bash

# New file to store the concatenated results
output_file="./mosdepth_summary.tsv"
> "$output_file" # Clear the output file if it already exists

# Global flag to check if the header has been written
header_written=0

# Iterate over each .bed file in the current directory
for file in *.bed; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract barcode identifier from the file name
        identifier=$(echo "$file" | sed -E 's/mosdepth\.(Barcode[0-9]+)\.regions\.bed/\1/')

        # Line counter
        line_count=0

        # Read the file line by line
        while IFS= read -r line; do
            # Increment line counter
            ((line_count++))

            # For the first file, include the header
            if [ $header_written -eq 0 ]; then
                if [ $line_count -eq 1 ]; then
                    # Append the new column name with a tab separator and set the header_written flag
                    echo -e "Chromosome\tStart\tEnd\tGene\tCoverage\tIdentifier" >> "$output_file"
                    header_written=1
                fi
            fi

            # Skip if line is empty
            [ -z "$line" ] && continue

            # For non-header lines, append the line with identifier information
            echo -e "$line\t$identifier" >> "$output_file"
        done < "$file"
    fi
done

echo "Concatenation complete. Output is in $output_file"
