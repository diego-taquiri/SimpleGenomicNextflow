#!/bin/bash

# New file to store the concatenated results
output_file="combined_summary.txt"
> "$output_file" # Clear the output file if it already exists

# Global flag to check if the header has been written
header_written=0

# Iterate over each *.kreport_bracken_species.txt in the current directory
for file in *.kreport_bracken_species.txt; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract barcode information from filename (assuming format barcodeXX.kreport_bracken_species.txt)
        barcode=$(echo "$file" | sed 's/.kreport_bracken_species.txt//')

        # Line counter
        line_count=0

        # Read the file line by line
        while IFS=$'\t' read -r percent covered_reads taxa_level rank id name; do
            # Increment line counter
            ((line_count++))

            # Skip the first line if it's a header or not needed
            if [ $line_count -eq 1 ]; then
                continue
            fi

            # Check if it's the first file's first line of data
            if [ $header_written -eq 0 ]; then
                # Append the new column name with a tab separator and set the header_written flag
                echo -e "Percent\tCovered Reads\tTaxa Level\tRank\tID\tName\tBarcode" >> "$output_file"
                header_written=1
            fi

            # For data lines, append the line with barcode information
            echo -e "$percent\t$covered_reads\t$taxa_level\t$rank\t$id\t$name\t$barcode" >> "$output_file"
        done < "$file"
    fi
done

echo "Concatenation complete. Output is in $output_file"
