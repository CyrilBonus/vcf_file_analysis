#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified on Jan 31, 2025
VCF Processing Script
"""
import pandas as pd
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

# File path to the VCF file
file = sys.argv[1]
# File path to the GenBank file
genbank_file = sys.argv[2]

def read_vcf(file):
    """
    Reads a VCF file and extracts relevant columns.
    """
    sections = []
    data = {}
    headers = []

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                if headers:
                    sections.append(pd.DataFrame(data))
                headers = line.strip().split('\t')
                data = {col: [] for col in headers}
                continue
            values = line.strip().split('\t')
            for col, val in zip(headers, values):
                data[col].append(val)
        if headers:
            sections.append(pd.DataFrame(data))

    return sections

def calculate_allele_frequencies(df):
    """
    Calculates allele frequencies from the INFO column.
    """
    allele_frequencies = []
    for info in df['INFO']:
        info_dict = {key: value for key, value in (item.split('=') for item in info.split(';') if '=' in item)}
        allele_frequencies.append(float(info_dict.get('AF', 'nan')))
    df['AF'] = allele_frequencies
    return df

def analyze_vcf(df, coding_regions):
    """
    Analyzes the VCF DataFrame and prints a summary.
    """
    total_variants = len(df)
    indels = df[df['INFO'].str.contains('SVTYPE=INS|SVTYPE=DEL')]
    high_quality_variants = df[df['QUAL'].astype(float) >= 30]
    passed_variants = df[df['FILTER'] == "PASS"]
    failed_variants = df[df['FILTER'] != "PASS"]

    # Calculate additional statistics
    avg_quality = df['QUAL'].astype(float).mean()
    quality_stats = df['QUAL'].astype(float).describe()
    df = calculate_allele_frequencies(df)
    allele_frequency_stats = df['AF'].describe()
    
    variants_per_chromosome = df['#CHROM'].value_counts()

    # Determine coding regions
    coding_variants = df[df.apply(lambda row: is_in_coding_region(int(row['POS']), coding_regions), axis=1)]
    coding_variants_count = len(coding_variants)
    coding_indels = coding_variants[coding_variants['INFO'].str.contains('SVTYPE=INS|SVTYPE=DEL')]

    return {
        "Total Variants": total_variants,
        "INDELs": len(indels),
        "High-Quality Variants": len(high_quality_variants),
        "Passed Variants": len(passed_variants),
        "Failed Variants": len(failed_variants),
        "Average Quality Score": avg_quality,
        "Quality Score Distribution": quality_stats.to_dict(),
        "Allele Frequency Distribution": allele_frequency_stats.to_dict(),
        "Variants per Chromosome": variants_per_chromosome.to_dict(),
        "Coding Variants": coding_variants_count,
        "Coding INDELs": len(coding_indels),
        "VCF Chromosomes": df['#CHROM'].unique().tolist(),
        "Positions": df['POS'].astype(int).tolist(),
        "Temperature": df['TEMPERATURE'].iloc[0] if 'TEMPERATURE' in df.columns else 'UNKNOWN'
    }

def get_coding_regions(genbank_file):
    """
    Parses the GenBank file to extract coding regions.
    """
    coding_regions = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start)
                end = int(feature.location.end)
                coding_regions.append((start, end))
    return coding_regions

def is_in_coding_region(pos, coding_regions):
    """
    Checks if a given position is within a coding region.
    """
    return any(start <= pos <= end for start, end in coding_regions)

def write_results(results):
    """
    Writes the analysis results to a file.
    """
    with open("results.txt", "w") as f:
        for section, result in results:
            f.write(f"{section} VCF File Analysis Results\n")
            f.write("=========================\n\n")
            f.write(f"Total Variants: {result['Total Variants']}\n")
            f.write(f"INDELs: {result['INDELs']}\n")
            f.write(f"High-Quality Variants (QUAL >= 30): {result['High-Quality Variants']}\n")
            f.write(f"Variants that Passed Filters: {result['Passed Variants']}\n")
            f.write(f"Variants that Failed Filters: {result['Failed Variants']}\n")
            f.write(f"Average Quality Score: {result['Average Quality Score']:.6f}\n\n")

            f.write("Quality Score Distribution:\n")
            for key, value in result['Quality Score Distribution'].items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")

            f.write("Allele Frequency Distribution:\n")
            for key, value in result['Allele Frequency Distribution'].items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")

            f.write("Variants per Chromosome:\n")
            for key, value in result["Variants per Chromosome"].items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")

            f.write(f"Coding Variants: {result['Coding Variants']}\n")
            f.write(f"Coding INDELs: {result['Coding INDELs']}\n")
            f.write("\n")

            f.write("Chromosome IDs in VCF file:\n")
            for chrom in result["VCF Chromosomes"]:
                f.write(f"  {chrom}\n")
            f.write("\n")

            # Generate and save the bar graph for variant positions
            generate_variant_position_graph(section, result['Positions'], result['Temperature'])

def generate_variant_position_graph(section, positions, temperature):
    """
    Generates a bar graph showing the positions of the variants in the chromosome.
    """
    plt.figure(figsize=(10, 6))
    plt.hist(positions, bins=50, edgecolor='black')
    plt.title(f'Variant Positions in Chromosome for {section} - Temperature: {temperature}')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{section}_variant_positions.png')
    plt.close()

# Example usage
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python vcf_file_analysis.py <vcf_file> <genbank_file>")
        sys.exit(1)

    sections = read_vcf(file)
    coding_regions = get_coding_regions(genbank_file)
    results = [(f"Section {i+1}", analyze_vcf(df, coding_regions)) for i, df in enumerate(sections)]
    write_results(results)