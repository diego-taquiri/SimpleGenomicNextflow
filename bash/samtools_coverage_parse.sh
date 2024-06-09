#!/bin/bash

# New file to store the concatenated results
output_file="concatenated_tables.txt"
> "$output_file" # Clear the output file if it already exists

# Global flag to check if the header has been written
header_written=0

# Iterate over files matching the pattern "barcode*_coverage.txt"
for file in barcode*_coverage.txt; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Get the filename without the extension
        filename=$(basename "$file" _coverage.txt)

        # Read the file line by line
        while IFS= read -r line; do
            # Skip if line is empty
            [ -z "$line" ] && continue

            # Check if it's the first file's header line
            if [ $header_written -eq 0 ]; then
                # Append the new column name with a tab separator and set the header_written flag
                echo -e "$line\tfilename" >> "$output_file"
                header_written=1
            else
                # For non-header lines or headers of subsequent files, append the line with filename
                # Skip the line if it's a header of a subsequent file
                [[ $line == *"#rname"* ]] && continue 
                echo -e "$line\t$filename" >> "$output_file"
            fi
        done < "$file"
    fi
done

echo "Concatenation complete. Output is in $output_file"
