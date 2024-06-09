#!/bin/bash

# Define base directory for NanoPlot results
NANOPLOT_DIR="./results/nanoplot"
SUMMARY_DIR="./results/summary"

# Ensure the summary directory exists
mkdir -p $SUMMARY_DIR

# Create an array with all samples based on the directories in the current path
samples=($(ls -d */ | cut -d'/' -f1))

# Define output table for NanoPlot summaries
SUMMARY_TABLE="$SUMMARY_DIR/nanoplot_summary_table.txt"

# Initialize the table with headers
echo -e "Sample\tMean Length\tMean Quality\tMedian Length\tMedian Quality\tNumber of Reads\tRead Length N50\tSTDEV Length\tTotal Bases" > $SUMMARY_TABLE

# Loop over each sample
for sample in "${samples[@]}"
do
    # Parse NanoStats.txt
    STAT_FILE="./$sample/NanoStats.txt"
    if [ -f "$STAT_FILE" ]; then
        echo "Parsing $STAT_FILE"
        MEAN_LENGTH=$(grep "Mean read length" $STAT_FILE | awk '{print $4}')
        MEAN_QUALITY=$(grep "Mean read quality" $STAT_FILE | awk '{print $4}')
        MEDIAN_LENGTH=$(grep "Median read length" $STAT_FILE | awk '{print $4}')
        MEDIAN_QUALITY=$(grep "Median read quality" $STAT_FILE | awk '{print $4}')
        NUMBER_READS=$(grep "Number of reads" $STAT_FILE | awk '{print $4}')
        READ_LENGTH_N50=$(grep "Read length N50" $STAT_FILE | awk '{print $4}')
        STDEV_LENGTH=$(grep "STDEV read length" $STAT_FILE | awk '{print $4}')
        TOTAL_BASES=$(grep "Total bases" $STAT_FILE | awk '{print $3}')
        
        if [ -n "$MEAN_LENGTH" ] && [ -n "$MEAN_QUALITY" ] && [ -n "$MEDIAN_LENGTH" ] && [ -n "$MEDIAN_QUALITY" ] && [ -n "$NUMBER_READS" ] && [ -n "$READ_LENGTH_N50" ] && [ -n "$STDEV_LENGTH" ] && [ -n "$TOTAL_BASES" ]; then
            echo -e "$sample\t$MEAN_LENGTH\t$MEAN_QUALITY\t$MEDIAN_LENGTH\t$MEDIAN_QUALITY\t$NUMBER_READS\t$READ_LENGTH_N50\t$STDEV_LENGTH\t$TOTAL_BASES" >> $SUMMARY_TABLE
        else
            echo "Failed to parse some values from $STAT_FILE" >&2
        fi
    else
        echo "NanoStats.txt not found for $sample" >&2
    fi
done
