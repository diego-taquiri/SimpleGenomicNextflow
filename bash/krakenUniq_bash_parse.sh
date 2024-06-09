#!/bin/bash

# New file to store the concatenated results
output_file="krakenUniq_summary.txt"
> "$output_file" # Clear the output file if it already exists

# Global flag to check if the header has been written
header_written=0

# Iterate over each report_file.tsv in barcode directories
for file in barcode*/report_file.tsv; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract barcode information from directory name
        barcode=$(dirname "$file")

        # Line counter
        line_count=0

        # Read the file line by line
        while IFS= read -r line; do
            # Increment line counter
            ((line_count++))

            # Skip the first two lines
            if [ $line_count -le 2 ]; then
                continue
            fi

            # Skip if line is empty
            [ -z "$line" ] && continue

            # Check if it's the first file's header line
            if [ $header_written -eq 0 ]; then
                # Append the new column name with a tab separator and set the header_written flag
                echo -e "$line\tbarcode" >> "$output_file"
                header_written=1
            else
                # For non-header lines or headers of subsequent files, append the line with barcode information
                # Skip the line if it's a header of a subsequent file
                [[ $line == *"Name"* ]] && continue 
                echo -e "$line\t$barcode" >> "$output_file"
            fi
        done < "$file"
    fi
done

echo "Concatenation complete. Output is in $output_file"
