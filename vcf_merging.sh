#!/bin/bash

# Output file name
output="combined.vcf"

# Check if output file already exists and remove it
if [ -f "$output" ]; then
    rm "$output"
fi

# Function to determine temperature based on filename
determine_temperature() {
    local filename=$1
    local sample=$(echo "$filename" | grep -oP 'P\d+-\d+')
    local number=$(echo "$sample" | grep -oP 'P\d+-(\d+)' | cut -d'-' -f2)
    if [ "$number" -ge 1 ] && [ "$number" -le 5 ]; then
        echo "COLD"
    else
        echo "HOT"
    fi
}

# Loop through all .vcf files in the current directory in natural order
for file in $(ls -v *.vcf); 
do
    # Determine the temperature for the current file
    temperature=$(determine_temperature "$file")
    filename_no_ext=$(echo "$file" | cut -d'.' -f1)
    echo "Processing file: $file, Temperature: $temperature"  # Debugging line
    
    # Add new header line before existing headers
    echo "##FILENAME=$filename_no_ext" >> "$output"
    
    # Process the file line by line
    while IFS= read -r line
    do
        if [[ $line == \#* ]]; then
            # If the line is a header line, add it to the output
            echo "$line" >> "$output"
        else
            # If the line is a data line, add the TEMPERATURE column
            echo -e "$line\t$temperature" >> "$output"
        fi
    done < "$file"
done

# Add the TEMPERATURE column name to the last header line
sed -i '/^#CHROM/s/$/\tTEMPERATURE/' "$output"

echo "All .vcf files have been merged into $output"
