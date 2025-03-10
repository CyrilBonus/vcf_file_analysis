#!/bin/bash

#ce script est utilisé que pour la fusion des vcf de la generation P15 car on ne les sépare pas en hot/cold 
#puisqu'ils n'ont pas encore subit de stress


if [ -z "$1" ]; then
    echo "ERREUR : pas de dossier spécifié donné"
    echo "il faut utiliser la commande suivante : ./fusion_vcf.sh /chemin/vers/dossier_contenant_les_vcf"
    exit 1
fi

data_dir="$1"  #dossier contenant les fichiers VCF
fusion_dir="${data_dir}/fusion_P15"  #dossier où stocker les fichiers fusionnés
mkdir -p "$fusion_dir"  #création du dossier fusion s'il n'existe pas

#fusion des fichiers P15
if [ -d "${data_dir}/P15" ]; then
    first_file=$(find "${data_dir}/P15" -name "P15-*.trimed1000.sv_sniffles.vcf" | head -n 1)

    if [ -n "$first_file" ]; then
        grep "^#" "$first_file" > "${fusion_dir}/P15.vcf"  #pour garder l'entete
        for file in "${data_dir}/P15"/*trimed1000.sv_sniffles.vcf; do
            if [ -f "$file" ]; then
                grep -v "^#" "$file" >> "${fusion_dir}/P15.vcf"
                echo "✅ Fusionné : $file"
            else
                echo "⚠️ Fichier manquant : $file"
            fi
        done
        echo "✨ fusion terminée pour P15.vcf ✨"
    else
        echo "aucun fichier trouvé dans ${data_dir}/P15/"
    fi
else
    echo "le répertoire ${data_dir}/P15 n'existe pas.."
fi

echo "✨ toutes les fusions sont terminées ! Les fichiers se trouvent dans : ${fusion_dir}/ "
