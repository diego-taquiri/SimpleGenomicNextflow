#!/bin/bash

# New file to store the concatenated results
output_file="checkm2_quality_summary.tsv"
> "$output_file" # Clear the output file if it already exists

# Global flag to check if the header has been written
header_written=0

# Iterate over each quality_report.tsv in each PER_UPCH* directory
for file in PER_UPCH*/quality_report.tsv; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract PER_UPCH directory name as an identifier
        identifier=$(dirname "$file" | cut -d'/' -f1)

        # Line counter
        line_count=0

        # Read the file line by line
        while IFS= read -r line; do
            # Increment line counter
            ((line_count++))

            # Skip the first line (headers in first file only)
            if [ $line_count -eq 1 ]; then
                if [ $header_written -eq 0 ]; then
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
