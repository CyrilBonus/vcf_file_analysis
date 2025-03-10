#!/usr/bin/env python

"""Ce script permet dans un premier temps d'extraire les informations des fichiers VCF et de les stocker dans un fichier CSV qui est plus facile d'utilisation.
Il permet aussi d'appliquer des filtres pour affiner l'analyse des variants en ne gardant que des variants pour lesquels on a une confiance élevée, 
pour cela les filtres utlisés sont les suivants : 
- QUAL > 30 : probabilité que l’appel du variant soit incorrect de 1/1000
- PASS tous les filtres 
- AF >= 0.1 : éviter le bruit
- SUPPORT > 10 : coverage minimum
 """

import os
import sys
import pandas as pd

#entrée et sortie
fusion_dir = sys.argv[1]   #dossier de fusion (issu du script 1)
output_dir = os.path.join(fusion_dir, "filtre")  
os.makedirs(output_dir, exist_ok=True)  #equivalent mkdir -p 

#filtrage et conversion en csv 
for file in os.listdir(fusion_dir):
    if file.endswith(".vcf"):  #que les fichiers VCF
        vcf_path = os.path.join(fusion_dir, file)
        output_csv = os.path.join(output_dir, file.replace(".vcf", "_filtre.csv"))

        #lire VCF + extraire données utiles
        filtered_variants = []
        with open(vcf_path, "r") as vcf:
            for line in vcf:
                if line.startswith("#"):  # Ignorer les commentaires
                    continue

                cols = line.strip().split("\t")
                chrom, pos, id_, ref, alt, qual, filt, info = cols[:8]
                info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in info.split(";") if "=" in kv}
                svlen = int(info_dict.get("SVLEN", 0))
                af = float(info_dict.get("AF", 0))
                support = int(info_dict.get("SUPPORT", 0))

                #filtres
                if float(qual) > 30 and filt == "PASS" and af >= 0.1 and support > 10:
                    #and svlen > 10 : je retire parce que ca m'enleve toutes les DEL sinon!
                    filtered_variants.append(cols[1:])  # On supprime seulement la colonne CHROM

        #convertir en CSV en gardant toutes les colonnes sauf CHROM
        df = pd.DataFrame(filtered_variants, columns=["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
        df.to_csv(output_csv, sep="\t", index=False)
        print(f"fichier filtré enregistré : {output_csv}")

print(f"tous les fichiers VCF fusionnés ont été filtrés : {output_dir}")
