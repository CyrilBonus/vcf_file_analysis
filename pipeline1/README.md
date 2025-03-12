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


## Code Functions

### 1. `read_vcf(file)`
**Description**:  
This function reads a VCF file, processes its content  and extracts what's considered relevant data columns for more analyses.

**Parameters**:  
- `file (str)`: The path to the VCF file.

**Returns**:  
- `list`: A list of tuples including the filename, the VCF DataFrame, and VCF type (meaning 'sniffles' or 'snp').

**Explanation**:  
This function opens and reads each line of a VCF file. It takes care of metadata lines, the column headers, and the data lines. For each section, it stores what's relevant in a pandas DataFrame, noticing the VCF type, and appends the section to a list. Finally, it returns a list of all sections in the file.

---

### 2. `detect_vcf_type(df)`
**Description**:  
This function determines whether the VCF file contains structural variants (Sniffles) or SNP-based variants.

**Parameters**:  
- `df (pd.DataFrame)`: The DataFrame containing VCF data.

**Returns**:  
- `str`: The type of VCF file (`'sniffles'` or `'snp'`).

**Explanation**:  
The function looks for the `SVTYPE` field in the `INFO` column. If it contains this field, it indicates the file contains structural variants (e.g., from Sniffles). If not, it assumes the file contains SNP-based variants.

---

### 3. `calculate_allele_frequencies(df)`
**Description**:  
This function calculates the allele frequencies from the `INFO` column in the VCF file and adds an 'AF' (Allele Frequency) column to the DataFrame.

**Parameters**:  
- `df (pd.DataFrame)`: The DataFrame containing VCF data.

**Returns**:  
- `pd.DataFrame`: The same DataFrame with an additional `AF` column representing allele frequencies.

**Explanation**:  
The function iterates over the `INFO` column, parses the allele frequencies (`AF`), and adds them as a new column in the DataFrame. If no allele frequency is found, it adds `NaN`.

---

### 4. `analyze_vcf(df, coding_regions)`
**Description**:  
This function analyzes a VCF DataFrame and generates summary statistics for the variants, including total counts, INDELs, quality scores, and allele frequencies.

**Parameters**:  
- `df (pd.DataFrame)`: The DataFrame containing VCF data.
- `coding_regions (list)`: A list of tuples representing coding regions from the GenBank file.

**Returns**:  
- `dict`: A dictionary with analysis results such as the total number of variants, INDEL counts, quality scores, allele frequencies, and more.

**Explanation**:  
This function performs several analyses:
- Counts the total number of variants.
- Filters variants based on type (INDELs), quality score, and filter status.
- Calculates the average quality score and provides a summary of the quality score distribution.
- Computes allele frequencies using the `calculate_allele_frequencies()` function.
- Determines the number of variants in coding regions by comparing positions with known coding regions from the GenBank file.

---

### 5. `get_coding_regions(genbank_file)`
**Description**:  
This function parses a GenBank file to extract coding regions (CDS) for comparison with VCF variant positions.

**Parameters**:  
- `genbank_file (str)`: The path to the GenBank file.

**Returns**:  
- `list`: A list of tuples containing the start and end positions of coding regions.

**Explanation**:  
This function reads the GenBank file, looks for features of type `CDS` (coding sequence), and extracts their start and end positions. It then returns a list of tuples containing these positions, which can be used to determine if a VCF variant falls within a coding region.

### 6. `is_in_coding_region(pos, coding_regions)`
**Description**:  
This function checks whether a given position falls within a coding region.

**Parameters**:  
- `pos (int)`: The position of the variant.
- `coding_regions (list)`: List of tuples containing coding region positions.

**Returns**:  
- `bool`: `True` if the position is within a coding region, `False` otherwise.

**Explanation**:  
The function iterates through the list of coding regions and checks if the given position (`pos`) falls within any of the start-end ranges. If so, it returns `True`; otherwise, `False`.

---

### 7. `generate_comparative_variant_graph_by_pnumber(results)`
**Description**:  
Generates a bar graph comparing the number of variants across all VCF files, grouped by 'PNUMBER'.

**Parameters**:  
- `results (list)`: A list of tuples containing the filename, VCF section, and analysis results.

**Explanation**:  
This function groups the results by 'PNUMBER' (e.g., 'P90', 'P15') and generates a bar graph showing the number of variants for each file in each group. It saves the graph as a PNG file for each PNUMBER group.

---

### 8. `calculate_mean_variants_by_pnumber(results)`
**Description**:  
Calculates the mean number of variants for each 'PNUMBER' group, considering both 'HOT' and 'COLD' files.

**Parameters**:  
- `results (list)`: A list of tuples containing the filename, VCF section, and analysis results.

**Returns**:  
- `dict`: A dictionary with the mean number of variants for each PNUMBER, separated by temperature.

**Explanation**:  
The function groups the results by 'PNUMBER' and temperature ('HOT' or 'COLD'), then calculates the mean number of variants for each group.

---

### 9. `generate_comparative_mean_variant_graph(mean_variants)`
**Description**:  
Generates a bar graph comparing the mean number of variants across different 'PNUMBER' groups, distinguishing between 'HOT' and 'COLD' files.

**Parameters**:  
- `mean_variants (dict)`: A dictionary with mean variant counts for each PNUMBER group.

**Explanation**:  
This function generates a bar graph comparing the mean variant counts for each PNUMBER, with separate bars for 'HOT' and 'COLD' temperatures. It saves the graph as a PNG file.

---

### 10. `write_results(results)`
**Description**:  
Writes the analysis results to a file.

**Parameters**:  
- `results (list)`: A list of tuples containing the filename, VCF section, and analysis results.

**Explanation**:  
This function writes detailed analysis results to a text file (`results.txt`), including variant counts, quality scores, allele frequencies, and more. It also generates additional graphs, such as variant position distributions and comparative graphs.

---

### 11. `generate_variant_position_graph(filename, section, positions, temperature)`
**Description**:  
Generates a bar graph showing the distribution of variant positions along the chromosome.

**Parameters**:  
- `filename (str)`: The name of the file used for the graph title and file name.
- `section (str)`: The section of the VCF file.
- `positions (list)`: A list of variant positions.
- `temperature (str)`: The temperature value from the analysis.

**Explanation**:  
This function creates a histogram showing the distribution of variant positions along the chromosome. It annotates the bars with the count values and saves the graph as a PNG file.  
