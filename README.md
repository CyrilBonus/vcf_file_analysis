# vcf_file_analysis

This project analyzes VCF files and generates several statistical outputs, that will be allele frequencies, variant counts, and the quality scores. This project uses Python3 and some bioinformatics libraries to process those VCF files and a reference genome.

## Scripts

1. **`vcf_merging.sh`**  
   Merges all `.vcf` files in the current directory into a single `combined.vcf` file.

2. **`vcf_cleanup.sh`**  
   Cleans up by removing the following files generated during the analysis:
   - `.txt` result files
   - `.png` graph files
   - The `combined.vcf` file

3. **`vcf_project.sh`**  
   Executes the main Python analysis code using the merged `combined.vcf` file and a reference GenBank file.  
   Command:  
   `./vcf_project.sh combined.vcf reference.gb`

## Requirements

- Python 3.x
- Pandas
- BioPython
- Matplotlib

## Usage

### Merging VCF Files
Run the `vcf_merging.sh` script to merge all the VCF files in the directory into a single `combined.vcf` file.

### Running the Analysis
Run the `vcf_project.sh` script, passing the merged `combined.vcf` file and a GenBank file as arguments.  
Example:  
`./vcf_project.sh combined.vcf reference.gb`

This will generate various outputs including graphs and a `results.txt` file containing summary statistics.

### Cleaning Up
Run the `vcf_cleanup.sh` script to remove all `.txt`, `.png`, and the `combined.vcf` file.
