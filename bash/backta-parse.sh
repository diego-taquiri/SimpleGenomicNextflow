#!/bin/bash

# New file to store the concatenated results
output_file="merged_assembly.tsv"
> "$output_file" # Clear the output file if it already exists

# Global flag to check if the header has been written
header_written=0

# Iterate over each assembly.tsv in each bakta.assembly.barcode* directory
for dir in bakta.assembly.barcode*/; do
    file="${dir}assembly.tsv"

    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract barcode directory name as an identifier
        identifier=$(basename "$dir")

        # Line counter
        line_count=0

        # Read the file line by line
        while IFS= read -r line; do
            # Increment line counter
            ((line_count++))

            # Skip the first 5 lines (headers in first file only)
            if [ $line_count -le 5 ]; then
                # For the first file, capture the header (line 6)
                if [ $line_count -eq 6 ] && [ $header_written -eq 0 ]; then
                    # Append the new column name with a tab separator and set the header_written flag
                    echo -e "$line\tIdentifier" >> "$output_file"
                    header_written=1
                fi
                continue
            fi

            # Skip if line is empty
            [ -z "$line" ] && continue

            # For non-header lines, append the line with identifier information
            echo -e "$line\t$identifier" >> "$output_file"
        done < "$file"
    fi
done

echo "Concatenation complete. Output is in $output_file"
