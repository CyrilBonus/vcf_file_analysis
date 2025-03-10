#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import statistics
import csv


def read_csv(file):
    """Parse the CSV file and extract variant information."""
    variants = []
    with open(file, 'r') as f:
        # First read the header line to get column names
        header = f.readline().strip().split()
        
        # Read the rest of the file
        for line in f:
            fields = line.strip().split()
            if len(fields) >= len(header):  # Ensure we have enough fields
                row = dict(zip(header, fields))
                variants.append({
                    "POS": int(row["POS"]),
                    "ID": row["ID"],
                    "REF": row["REF"],
                    "ALT": row["ALT"],
                    "QUAL": row["QUAL"],
                    "FILTER": row["FILTER"],
                    "INFO": row["INFO"],
                    "FORMAT": row["FORMAT"],
                    "SAMPLE": row["SAMPLE"]
                })
    return variants


def count_variants_by_type(variants, variant_types=["INS", "DEL", "DUP"]):
    """Count variants based on SVTYPE in INFO field."""
    counts = {var_type: 0 for var_type in variant_types}
    for variant in variants:
        for var_type in variant_types:
            if f"SVTYPE={var_type}" in variant["INFO"]:
                counts[var_type] += 1
    return counts


def count_total_variants(variants):
    """Count total variants."""
    return len(variants)


def count_variants_by_filter(variants):
    """Count variants by FILTER value."""
    filter_counts = {}
    for variant in variants:
        filter_counts[variant["FILTER"]] = filter_counts.get(variant["FILTER"], 0) + 1
    return filter_counts


def allele_frequency_distribution(variants):
    """Calculate statistics of allele frequencies (AF) from INFO fields."""
    af_values = []
    for variant in variants:
        for field in variant["INFO"].split(";"):
            if field.startswith("AF="):
                try:
                    af_values.append(float(field.split("=")[1]))
                except ValueError:
                    pass
    if af_values:
        mean_af = statistics.mean(af_values)
        stdev_af = statistics.stdev(af_values) if len(af_values) > 1 else 0
    else:
        mean_af, stdev_af = 0, 0

    return mean_af, stdev_af, len(af_values)


def write_summary_to_csv(output_path, counts_by_type, total_variants, filter_counts, af_stats):
    """Write analysis summary to a CSV file."""
    mean_af, stdev_af, af_count = af_stats
    with open(output_path, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Total Variants", total_variants])
        
        writer.writerow(["\nVariant Type Counts", ""])
        for var_type, count in counts_by_type.items():
            writer.writerow([var_type, count])

        writer.writerow(["\nFilter Counts", ""])
        for filt, count in filter_counts.items():
            writer.writerow([filt, count])

        writer.writerow(["\nAllele Frequency Distribution", ""])
        writer.writerow(["Mean AF", f"{mean_af:.3f}"])
        writer.writerow(["Standard Deviation AF", f"{stdev_af:.3f}"])
        writer.writerow(["AF Count", af_count])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <directory_path>")
        sys.exit(1)

    directory = sys.argv[1]

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    # Process all CSV files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".csv") and not filename.endswith("_stats.csv"):
            file_path = os.path.join(directory, filename)
            try:
                variants = read_csv(file_path)

                # Compute metrics
                counts_by_type = count_variants_by_type(variants)
                total_variants = count_total_variants(variants)
                filter_counts = count_variants_by_filter(variants)
                af_stats = allele_frequency_distribution(variants)

                # Generate output file path
                output_file = os.path.join(directory, f"{os.path.splitext(filename)[0]}_stats.csv")
                write_summary_to_csv(output_file, counts_by_type, total_variants, filter_counts, af_stats)

                print(f"Analysis for {filename} completed: {output_file}")
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")

    print("All CSV files processed.")