#!/bin/bash

vcf_file="$1"
reference_file="$2"

# part 1: check if the files are actual VCF and reference genome files before running the rest

controle=false

if [ -f "$vcf_file" ]; then
    if [[ "$vcf_file" == *.vcf ]]; then
        echo "'$vcf_file' is a VCF file"
        controle=true
    else
        echo "'$vcf_file' is a regular file but not a VCF file"
    fi
else
    echo "'$vcf_file' is not a regular file or does not exist"
fi

if [ -f "$reference_file" ]; then
    echo "'$reference_file' is a reference genome file"
else
    echo "'$reference_file' is not a regular file or does not exist"
    controle=false
fi

# part 2: if the files pass the control, run the python analysis

if [ "$controle" = true ]; then
    python3 ./VCF_file_analysis.py "$vcf_file" "$reference_file" > results.txt
else
    echo "One or both files are not valid; we cannot run the analysis"
fi
