#!/bin/bash

# Loop through files matching the pattern *_N*_2*
for l in *_N*_2*; do
    echo "Processing file: $l"

    # Extract the folder name from the first 6 characters of the file name
    foldername=$(echo "$l" | cut -c1-6)

    # Create a new directory named after the extracted folder name
    mkdir -p "$foldername"

    # Move the file into the newly created directory
    mv "$l" "$foldername"

    # Change to the new directory
    cd "$foldername" || exit

    # Run cutadapt to trim adapters, discarding untrimmed sequences
    cutadapt -a GATGTCCACGAGGTCTCT...CGTACGCTGCAGGTCGAC -o trimmed_"$foldername"_r2.fastq.gz --discard-untrimmed "$l"

    # Keep only the sequence of the barcode
    zcat trimmed_"$foldername"_r2.fastq.gz | paste - - - - | cut -f2 > "$foldername"_r2_seq.txt

    # Count the number of barcodes and store in a variable
    reads=$(wc -l < "$foldername"_r2_seq.txt)
    echo "Number of barcodes for $foldername: $reads"

    # Optionally, add this to a summary file if needed
    echo "$foldername: $reads" >> ../barcode_summary.txt

    # Calculate read lengths and save to read_length_<foldername>.txt
    zcat trimmed_"$foldername"_r2.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > "read_length_${foldername}.txt"

    # Return to the parent directory
    cd ..
done
