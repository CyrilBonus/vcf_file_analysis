import pandas as pd
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

# File path to the VCF file provided as a command-line argument
file = sys.argv[1]
# File path to the GenBank file provided as a command-line argument
genbank_file = sys.argv[2]

def read_vcf(file):
    """
    Reads a VCF (Variant Call Format) file, processes its content, 
    and extracts relevant data columns for further analysis.
    
    Args:
        file (str): Path to the VCF file.
        
    Returns:
        list: A list of tuples containing filename, VCF DataFrame, and VCF type ('sniffles' or 'snp').
    """
    sections = []
    data = {}
    headers = []
    filename = None  # Variable to track the current filename in the VCF

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('##'):  # Skip metadata lines
                if line.startswith('##FILENAME='):
                    if filename and headers:  # Save data from previous section
                        try:
                            df = pd.DataFrame(data)
                            vcf_type = detect_vcf_type(df)
                            sections.append((filename, df, vcf_type))
                        except ValueError as e:
                            print(f"Skipping section {filename} due to inconsistent data: {e}")
                    filename = line.strip().split('=')[1]  # Extract the filename from metadata
                    data = {}  # Reset data dictionary for the new file
                    headers = []  # Reset headers for the new file
                continue

            if line.startswith('#'):  # Header line (column names)
                if headers and data:  # Ensure the previous section is stored before starting a new one
                    try:
                        df = pd.DataFrame(data)
                        vcf_type = detect_vcf_type(df)
                        sections.append((filename, df, vcf_type))
                    except ValueError as e:
                        print(f"Skipping section {filename} due to inconsistent data: {e}")
                headers = line.strip().split('\t')  # Parse the column headers
                data = {col: [] for col in headers}  # Prepare an empty data dictionary for new columns
                continue

            values = line.strip().split('\t')  # Data line (VCF entries)
            if len(values) == len(headers):  # Ensure correct number of columns
                for col, val in zip(headers, values):  # Add each value to the appropriate column
                    data[col].append(val)
            else:
                print(f"Skipping line due to inconsistent column length: {line.strip()}")

        # Ensure the last section is not skipped
        if headers and data:
            try:
                df = pd.DataFrame(data)
                vcf_type = detect_vcf_type(df)
                sections.append((filename, df, vcf_type))
            except ValueError as e:
                print(f"Skipping section {filename} due to inconsistent data: {e}")

    return sections

def detect_vcf_type(df):
    """
    Determines whether the VCF file is from Sniffles (structural variants, SVTYPE) or SNP-based.
    
    Args:
        df (pd.DataFrame): DataFrame containing VCF data.
        
    Returns:
        str: The type of VCF file ('sniffles' or 'snp').
    """
    if df['INFO'].str.contains('SVTYPE=').any():
        return "sniffles"  # Structural variants (Sniffles)
    return "snp"  # SNP-based file

def calculate_allele_frequencies(df):
    """
    Calculates the allele frequencies from the 'INFO' column in the VCF file.
    
    Args:
        df (pd.DataFrame): DataFrame containing VCF data.
        
    Returns:
        pd.DataFrame: DataFrame with added 'AF' column representing allele frequencies.
    """
    allele_frequencies = []
    for info in df['INFO']:
        info_dict = {key: value for key, value in (item.split('=') for item in info.split(';') if '=' in item)}
        allele_frequencies.append(float(info_dict.get('AF', 'nan')))  # Retrieve AF or NaN if not present
    df['AF'] = allele_frequencies  # Add the allele frequencies to the DataFrame
    return df

def analyze_vcf(df, coding_regions):
    """
    Analyzes the VCF DataFrame and generates summary statistics for the variants.
    
    Args:
        df (pd.DataFrame): DataFrame containing VCF data.
        coding_regions (list): List of tuples representing coding regions from the GenBank file.
        
    Returns:
        dict: A dictionary with analysis results such as the total number of variants, INDEL counts, quality scores, etc.
    """
    total_variants = len(df)
    indels = df[df['INFO'].str.contains('SVTYPE=INS|SVTYPE=DEL')]  # Filter for insertion/deletion variants
    high_quality_variants = df[df['QUAL'].astype(float) >= 30]  # Filter for variants with a quality score >= 30
    passed_variants = df[df['FILTER'] == "PASS"]  # Filter for variants that passed filters
    failed_variants = df[df['FILTER'] != "PASS"]  # Filter for variants that failed filters

    # Calculate additional statistics
    avg_quality = df['QUAL'].astype(float).mean()  # Average quality score
    quality_stats = df['QUAL'].astype(float).describe()  # Summary of quality score distribution
    df = calculate_allele_frequencies(df)  # Calculate allele frequencies
    allele_frequency_stats = df['AF'].describe()  # Summary of allele frequency distribution
    
    variants_per_chromosome = df['#CHROM'].value_counts()  # Count variants per chromosome

    # Determine coding regions by checking if variant positions overlap with known coding regions
    coding_variants = df[df.apply(lambda row: is_in_coding_region(int(row['POS']), coding_regions), axis=1)]
    coding_variants_count = len(coding_variants)
    
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
        "VCF Chromosomes": df['#CHROM'].unique().tolist(),
        "Positions": df['POS'].astype(int).tolist(),
        "Temperature": df['TEMPERATURE'].iloc[0] if 'TEMPERATURE' in df.columns else 'UNKNOWN'
    }

def get_coding_regions(genbank_file):
    """
    Parses the GenBank file to extract coding regions (CDS) for comparison with VCF variant positions.
    
    Args:
        genbank_file (str): Path to the GenBank file.
        
    Returns:
        list: A list of tuples containing the start and end positions of coding regions.
    """
    coding_regions = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":  # Extract coding sequence features
                start = int(feature.location.start)
                end = int(feature.location.end)
                coding_regions.append((start, end))  # Add each coding region to the list
    return coding_regions

def is_in_coding_region(pos, coding_regions):
    """
    Checks if a given position is within a coding region.
    
    Args:
        pos (int): The position of the variant.
        coding_regions (list): List of tuples containing coding region positions.
        
    Returns:
        bool: True if the position is within a coding region, False otherwise.
    """
    return any(start <= pos <= end for start, end in coding_regions)

import re

def generate_comparative_variant_graph_by_pnumber(results):
    """
    Generates a bar graph comparing the number of variants across all VCF files, grouped by 'PNUMBER'.
    
    Args:
        results (list): A list of tuples containing the filename, VCF section, and analysis results.
    """
    # Group files by 'PNUMBER' (e.g., 'P90', 'P15')
    grouped_results = {}
    for filename, section, result in results:
        pnumber_match = re.match(r'^(P\d{2})-', filename)  # Match the PNUMBER prefix (e.g., 'P90')
        if pnumber_match:
            pnumber = pnumber_match.group(1)  # Extract the PNUMBER
            if pnumber not in grouped_results:
                grouped_results[pnumber] = []
            grouped_results[pnumber].append((filename, result))  # Group results by PNUMBER
    
    # Generate a graph for each group
    for pnumber, group in grouped_results.items():
        filenames = [filename for filename, _ in group]
        total_variants = [result["Total Variants"] for _, result in group]

        plt.figure(figsize=(10, 6))
        plt.bar(filenames, total_variants, color='purple', edgecolor='black')
        plt.xlabel("VCF Files")
        plt.ylabel("Total Variants")
        plt.title(f"Comparison of Variant Counts for {pnumber}")
        plt.xticks(rotation=45, ha="right")
        plt.grid(axis="y", linestyle="--", alpha=0.7)

        # Annotate bars with values
        for i, count in enumerate(total_variants):
            plt.text(i, count + 2, str(count), ha="center", fontsize=10)

        plt.tight_layout()
        plt.savefig(f"variant_comparison_{pnumber}.png")
        plt.close()

def calculate_mean_variants_by_pnumber(results):
    """
    Calculates the mean number of variants for each 'PNUMBER' group.
    
    Args:
        results (list): A list of tuples containing the filename, VCF section, and analysis results.
        
    Returns:
        dict: A dictionary with the mean number of variants for each PNUMBER.
    """
    # Group results by 'PNUMBER'
    grouped_results = {}
    for filename, section, result in results:
        pnumber_match = re.match(r'^(P\d{2})-', filename)
        if pnumber_match:
            pnumber = pnumber_match.group(1)
            if pnumber not in grouped_results:
                grouped_results[pnumber] = []
            grouped_results[pnumber].append(result["Total Variants"])
    
    # Calculate mean for each PNUMBER group
    mean_variants = {pnumber: sum(variants) / len(variants) for pnumber, variants in grouped_results.items()}
    return mean_variants

def generate_comparative_mean_variant_graph(mean_variants):
    """
    Generates a bar graph comparing the mean number of variants across different 'PNUMBER' groups.
    
    Args:
        mean_variants (dict): A dictionary with mean variant counts for each PNUMBER group.
    """
    pnumbers = list(mean_variants.keys())
    means = list(mean_variants.values())

    plt.figure(figsize=(10, 6))
    plt.bar(pnumbers, means, color='green', edgecolor='black')
    plt.xlabel("PNUMBER")
    plt.ylabel("Mean Number of Variants")
    plt.title("Comparison of Mean Variant Counts Across PNUMBER Groups")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Annotate bars with values
    for i, mean in enumerate(means):
        plt.text(i, mean + 2, f'{mean:.2f}', ha="center", fontsize=10)

    plt.tight_layout()
    plt.savefig("mean_variant_comparison.png")
    plt.close()

def write_results(results):
    """
    Writes the analysis results to a file.
    """
    with open("results.txt", "w") as f:
        for filename, section, result in results:
            # Debugging: Print to ensure results are being processed

            f.write(f"{filename} - {section} VCF File Analysis Results\n")
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

            f.write(f"Coding Variants: {result['Coding Variants']}/{result['Total Variants']}\n")
            f.write("\n")

            f.write("Chromosome IDs in VCF file:\n")
            for chrom in result["VCF Chromosomes"]:
                f.write(f"  {chrom}\n")
            f.write("\n")

            # Generate and save the bar graph for variant positions
            generate_variant_position_graph(filename, section, result['Positions'], result['Temperature'])
            
        # Calculate the mean number of variants per PNUMBER and generate the comparative graph
        mean_variants = calculate_mean_variants_by_pnumber(results)
        generate_comparative_mean_variant_graph(mean_variants)
        generate_comparative_variant_graph_by_pnumber(results)




def generate_variant_position_graph(filename, section, positions, temperature):
    """
    Generates a bar graph showing the distribution of variant positions along the chromosome.
    
    Args:
        filename (str): The name of the file to be used in the graph title and file name.
        section (str): The section of the VCF file.
        positions (list): A list of variant positions.
        temperature (str): The temperature value from the analysis.
    """
    plt.figure(figsize=(10, 6))
    # Generate the histogram and store counts, bins, and patches
    counts, bins, patches = plt.hist(positions, bins=50, edgecolor='black', color='blue')

    # Annotate max values on top of each bar
    for count, bin_edge in zip(counts, bins[:-1]):  
        if count > 0:  # Only label non-zero bars
            plt.text(bin_edge + (bins[1] - bins[0]) / 2, count, f'{int(count)}', ha='center', va='bottom', fontsize=9, rotation=45)

    # Set plot title and labels
    plt.title(f'Variant Positions in Chromosome for {filename} - {section} - Temperature: {temperature}')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.tight_layout()

    # Save the figure as a PNG file with a unique name
    plt.savefig(f'{filename}_{section}_variant_positions.png')
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python vcf_file_analysis.py <vcf_file> <genbank_file>")
        sys.exit(1)

    # Assign file paths from command line arguments
    file = sys.argv[1]
    genbank_file = sys.argv[2]

    # Read VCF file and extract data sections
    sections = read_vcf(file)
    
    # Parse GenBank file for coding regions
    coding_regions = get_coding_regions(genbank_file)
    
    results = []

    # Analyze each VCF section and store results
    for filename, df, vcf_type in sections:
        section_name = f"Section {vcf_type}"
        result = analyze_vcf(df, coding_regions)
        results.append((filename, section_name, result))

    # Write results to a file
    write_results(results)
