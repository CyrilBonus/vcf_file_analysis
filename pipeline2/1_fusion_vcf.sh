#!/bin/bash



#vérification du dossier donné en argument
if [ -z "$1" ]; then
    echo "ERREUR : pas de dossier spécifié donné"
    echo "il faut utiliser la commande suivante : ./fusion_vcf.sh /chemin/vers/dossier_contenant_les_vcf"
    exit 1
fi

data_dir="$1"  #dossier contenant les fichiers VCF
fusion_dir="${data_dir}/fusion"  #dossier où stocker les fichiers fusionnés
mkdir -p "$fusion_dir"  #création du dossier fusion s'il n'existe pas

#fusion pour chaque passage
for passage in P65; do
    echo "fusion des fichiers pour $passage..."

    #fichiers froids (1 à 5)
    first_file=$(ls ${data_dir}/${passage}/${passage}-1.trimed1000.sv_sniffles.vcf | head -n 1)
    if [ -f "$first_file" ]; then
        grep "^#" "$first_file" > ${fusion_dir}/${passage}_cold.vcf  #garder l'entête
        cat ${data_dir}/${passage}/${passage}-{1..5}.trimed1000.sv_sniffles.vcf | grep -v "^#" >> ${fusion_dir}/${passage}_cold.vcf
        echo "fusion terminée pour ${passage}_cold.vcf"
    else
        echo "pas de fichier froid trouvé pour ${passage}, fusion ignorée."
    fi

    #fichiers chauds (6 à 10)
    first_file=$(ls ${data_dir}/${passage}/${passage}-6.trimed1000.sv_sniffles.vcf | head -n 1)
    if [ -f "$first_file" ]; then
        grep "^#" "$first_file" > ${fusion_dir}/${passage}_hot.vcf  #garder l'entête
        cat ${data_dir}/${passage}/${passage}-{6..10}.trimed1000.sv_sniffles.vcf | grep -v "^#" >> ${fusion_dir}/${passage}_hot.vcf
        echo "fusion terminée pour ${passage}_hot.vcf"
    else
        echo "pas de fichier chaud trouvé pour ${passage}, fusion ignorée."
    fi
done

echo ""✨ toutes les fusions sont terminées ! Les fichiers se trouvent dans : ${fusion_dir}/ "


