import re
import sys
from collections import defaultdict
import matplotlib.pyplot as plt

def count_variants(vcf_file, output_file):
    variant_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # PNUMBER -> sample -> temperature -> count
    
    with open(vcf_file, 'r') as file, open(output_file, 'w') as output:
        current_sample = None
        current_pnumber = None
        
        for line in file:
            line = line.strip()
            
            # Extract PNUMBER and sample ID from filename header
            if line.startswith("##FILENAME="):
                match = re.search(r'##FILENAME=([^\.]+)', line)
                if match:
                    full_sample_id = match.group(1)
                    base_sample = full_sample_id.split('-')[0]  # Extract PNUMBER (e.g., P90 from P90-1)
                    current_sample = full_sample_id
                    
                    # Track the current PNUMBER to sort later
                    if current_pnumber != base_sample:
                        current_pnumber = base_sample
                    # Don't write debug prints, so skipping here
                    
            # Process variant lines
            elif not line.startswith("#") and current_sample:
                fields = line.split("\t")
                if len(fields) < 9:
                    continue  # Skip malformed lines
                
                temperature = fields[-1]  # Last field is temperature
                
                # Only count if temperature matches valid ones
                if temperature in ['COLD', 'HOT', 'NORMAL']:
                    base_sample = current_sample.split('-')[0]  # Extract PNUMBER
                    
                    # Here, we check if the variant is present by looking at the GT field (0/1, 1/1, etc.)
                    genotypes = fields[9:]  # Genotypes are in columns after the 9th field
                    
                    variant_found = False
                    for genotype in genotypes:
                        if '1' in genotype:  # If there is an alternate allele, it's a variant
                            variant_found = True
                            break
                    
                    if variant_found:
                        variant_counts[base_sample][current_sample][temperature] += 1
    
        # Write the variant counts to the output file, ensuring the samples are sorted numerically
        for p_number, samples in sorted(variant_counts.items(), key=lambda x: int(x[0].split('-')[0][1:])):
            output.write(f"\n{p_number} Variant Counts:\n")
            # Sort the samples numerically
            for sample, temps in sorted(samples.items(), key=lambda x: int(x[0].split('-')[1])):
                output.write(f"  Sample: {sample}\n")
                for temp, count in temps.items():
                    output.write(f"    {temp}: {count} variants\n")
    
    return variant_counts

def plot_variants(variant_counts):
    for p_number, samples in variant_counts.items():
        plt.figure(figsize=(10, 6))

        # Sort sample names numerically
        sample_names = sorted(samples.keys(), key=lambda x: int(x.split('-')[1]))

        # Debug: Check the sorting order of the sample names
        print(f"Sorted sample names for {p_number}: {sample_names}")

        temperatures = ['COLD', 'HOT', 'NORMAL']
        
        # Prepare data for plotting with the correctly sorted sample names
        data = {temp: [samples[sample].get(temp, 0) for sample in sample_names] for temp in temperatures}
        
        x = range(len(sample_names))
        bar_width = 0.25
        
        # Create bars for each temperature
        for i, temp in enumerate(temperatures):
            plt.bar(x, data[temp], width=bar_width, label=temp, align='center')  # Use 'align=center'
        
        plt.xlabel('Samples')
        plt.ylabel('Number of Variants')
        plt.title(f'Variants per Sample in {p_number}')
        
        # Set x-ticks based on sorted sample names
        plt.xticks(x, sample_names, rotation=45, ha='right')
        
        plt.legend()
        plt.tight_layout()

        # Save the plot and close
        plt.savefig(f'{p_number}_variants.png')
        plt.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <vcf_file> <output_file>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    output_file = sys.argv[2]  # Take output file name as an argument
    variant_counts = count_variants(vcf_file, output_file)
    plot_variants(variant_counts)
