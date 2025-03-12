import os
import sys
import pandas as pd
from collections import defaultdict

#objectif de ce script : cat les csv
#- ajouter la provenance de chaque variant grace au nom du fichier dont il est extrait dans une nouvelle colonne
#- garder chaque ligne distincte : c'est a dire qu'il y a peut etre plusieurs fois le meme variant dans le meme fichier, je veux pas

#plus tard je veux changer la maniere dont la provenance est spécifiée : 
#- transformer le string "P15_P90" en deux cases 
    #- passage start : 15
    #- passage end : 90

#cat les variants identiques a partir de variant_ID = position, svtype, svlen... mais pour passage start et end : garder start min et end max


#lecture 
def process_csv_files(directory):
    variant_dict = defaultdict(list)
    all_data = []

    #boucle sur les fichiers csv du repertoire
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            filepath = os.path.join(directory, filename)
            try:
                #vide?
                if os.path.getsize(filepath) == 0:
                    print(f"Warning: File {filename} is empty. Skipping.")
                    continue
                df = pd.read_csv(filepath)
                df_ORF = df[df['ORF'].notna() & (df['ORF'] != "") & (df['ORF'] != "None")] #que les mutations sur ORF

                if not df_ORF.empty:
                    all_data.append(df_ORF)

            except pd.errors.EmptyDataError:
                print(f"Warning: File {filename} is empty or has no valid data. Skipping.")
            except Exception as e:
                print(f"Error processing file {filename}: {e}. Skipping.")

    # Combine all data into a single DataFrame
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)

        return combined_df
    else:
        print("No valid data found in any files.")
        return pd.DataFrame()

# Directory containing the CSV files
directory = sys.argv[1]

# Process the CSV files and get the combined DataFrame
combined_df = process_csv_files(directory)

# Save the results to a new CSV file
if not combined_df.empty:
    output_file = f"{directory}/all_combined.csv"
    combined_df.to_csv(output_file, index=False)
    print(f" ✨Yay!! We combined the variants from {directory} into {output_file}✨")
else:
    print("No data to save.")