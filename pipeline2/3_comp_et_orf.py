#!/usr/bin/env python

"""
Ce script permet dans un premier temps de comparer les variants présents entre différents passage pour déterminer lesquels sont conservés, disparaissent ou apparaissent entre deux générations différentes.
Ensuite, grace a une lecture du fichier fasta des ORF de la référence il permet de trouver les variations touchant des ORF annotées.
Le script a été conçu et réalisé en collaboration avec Hermine Kioussou et utilise les fonctions de comparaison de Elio Torquet.
Ce script n'est pas encore automatisé, il est nécessaire de remplacer les générations à chaque utilisation en veillant à respecter l'ordre chronologique des génération. 
Cette version permet la comparaison entre P50 et P65, pour une comparaison des passages 50 et 65, veiller à bien remplacer P50 par P50 et P65 par P65 (partout dans le script)!
 """

import os
import sys
import pandas as pd
from Bio import SeqIO
import re

##Partie 1 : initialisation 
print("1 - Initialisation du script...")

#définition du dossier où sont stockés les fichiers filtrés
if len(sys.argv) > 1:
    dossier_filtrage = sys.argv[1]
else:
    dossier_filtrage = "./filtrage"  # Valeur par défaut si aucun argument n'est donné
if not os.path.exists(dossier_filtrage):
    print(f"ERREUR : le dossier {dossier_filtrage} n'existe pas !")
    exit()

##Partie 2 : lecture et stockage des variants
print(f"2 - Lecture des fichiers dans le dossier : {dossier_filtrage}")
fichiers = [f for f in os.listdir(dossier_filtrage) if f.endswith(".csv")] #os.listdir(dossier_filtrage) permet de lister tous les fichiers présents dans le dossier spécifié (dossier_filtrage)
if not fichiers:
    print("pas de fichier CSV trouvé dans le dossier filtrage !")
    exit()

mutations = {}
conditions = ["hot", "cold"]

#stocker les variants dans un df
for fichier in fichiers:
    chemin_fichier = os.path.join(dossier_filtrage, fichier)
    print(f"Lecture de {fichier}...")
    try:
        df = pd.read_csv(chemin_fichier, sep="\t") 
        print(f"Colonnes disponibles dans {fichier} : {df.columns.tolist()}")
    except Exception as e:
        print(f" Erreur de lecture pour {fichier} : {e}")
        continue
    #vérifier la colonne INFO avant d'extraire les valeurs
    if "INFO" not in df.columns:
        print(f" Fichier {fichier} ignoré : colonne INFO manquante !")
        continue
    #si une colonne essentiel info (qui contient ["pos", "len", "svtype", "svlen", "alt", "ref", "coverage", "Af"]) 
    #manque on affiche un message d'erreur et on continue
    
    #détermination de la condition et du passage
    condition = "hot" if "hot" in fichier else "cold" 
    passage = "P50" if "P50" in fichier else "P65"
    
    if passage not in mutations:
        mutations[passage] = {}
    if condition not in mutations[passage]: #c'est ici on fait le stockage dans mutation 
        mutations[passage][condition] = []
    
    for _, row in df.iterrows(): #parcourt chaque ligne du fichier CSV etstockage de chaque mutation sous forme de dictionnaire
        try : 
            info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in row["INFO"].split(";") if "=" in kv} #séparation des éléments de la colonne info       
            variant = {
                "pos": int(row["POS"]), 
                "svtype": info_dict.get("SVTYPE", "UNKNOWN"),
                "svlen": int(info_dict.get("SVLEN", 0)),
                "alt": row["ALT"],
                "ref": row["REF"],
                "coverage": int(info_dict.get("COVERAGE", "0").split(",")[0]),
                "Af": float(info_dict.get("AF", "0.0").split(",")[0])
            }
            mutations[passage][condition].append(variant)

        except Exception as e:
            print(f"ERREUR lors de l'extraction des données pour une ligne dans {fichier} : {e}")
            continue    
print("Lecture des fichiers terminée !")

##Partie 3 : comparaison des variants entre les deux passages

#calculer l'identité de séquence
def seq_identity(s1, s2):
    identity = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            identity += 1
    return identity / len(s1)

#comparer deux variants
def variant_equal(v1, v2, sim_thresold=1):
    if v1["svtype"] != v2["svtype"]: 
        #on compare deux mutations pour voir si elles sont du même type : 
        #les délétions avec les délétions, les insertions avec les insertions
        return False
    
    full_length = max(v1["pos"] + abs(v1["svlen"]), v2["pos"] + abs(v2["svlen"])) - min(v1["pos"], v2["pos"])
    common_length = max(0, min(v1["pos"] + abs(v1["svlen"]), v2["pos"] + abs(v2["svlen"])) - max(v1["pos"], v2["pos"]))
    shared = common_length / full_length if full_length > 0 else 0

    if v1["svtype"] == "INS" and shared > 0 and v1["alt"] != "<INS>" and v2["alt"] != "<INS>":
        first, second = (v1, v2) if v1["pos"] < v2["pos"] else (v2, v1)
        common_start = second["pos"] - first["pos"]
        seq1, seq2 = first["alt"][common_start:], second["alt"]
        common_stop = min(len(seq1), len(seq2))
        shared *= seq_identity(seq1[:common_stop], seq2[:common_stop])
    
    return shared >= sim_thresold

print("3 - Comparaison des mutations entre P50 et P65...")

mutations_apparues = [] #liste pour stocker les mutations suivant le status 
mutations_disparues = []
mutations_conservees = []

for condition in conditions: #conditions défini en ligne 34
    if "P50" not in mutations or "P65" not in mutations:
        print("ERREUR : Les mutations pour P50 ou P65 sont manquantes !") 
        #Si P50 ou P65 n'existe pas dans le dictionnaire mutations, 
        #alors on ne peut pas faire la comparaison et on arrête le script (exit()).
        exit()
    
    mutations_P50 = mutations["P50"].get(condition, [])
    mutations_P65 = mutations["P65"].get(condition, [])
    
    for variant_P65 in mutations_P65:
        found = any(variant_equal(variant_P65, v) for v in mutations_P50) 
        #on cherche si le variant existe déjà dans P50 grâce a la fonction variant_equal()
        if found:
            mutations_conservees.append({"condition": condition, **variant_P65}) 
            #si oui on l'ajoute à la liste mutation_conservees
        else:
            mutations_apparues.append({"condition": condition, **variant_P65}) 
            #si non la mutation est apparue à P65
    
    for variant_P50 in mutations_P50:
        found = any(variant_equal(variant_P50, v) for v in mutations_P65) 
        #on cherche si la mutation de P50 existe encore à P65
        if not found:
            mutations_disparues.append({"condition": condition, **variant_P50}) 
            #si non on l'ajoute à la liste mutations_disparues 
print("Comparaison terminée !")


##Partie 4 : lecture du fichier fasta de la reference pour stocker les ORF

orf_data = []
fichier_orf = input("Entrez le chemin du fichier ORF.fasta s'il vous plait : ")
if os.path.isdir(fichier_orf):  #vérifie si c'est un dossier
    fichiers_fasta = [f for f in os.listdir(fichier_orf) if f.endswith(".fasta")]
    if fichiers_fasta:
        fichier_orf = os.path.join(fichier_orf, fichiers_fasta[0])  #prend le premier fichier trouvé
        print(f"fichier ORF trouvé : {fichier_orf}")
    else:
        print("ERREUR : Aucun fichier .fasta trouvé dans le dossier.")
        exit()
if not os.path.exists(fichier_orf): #verifier si le fichier existe
    print("ERREUR : Le fichier ORF spécifié n'existe pas.")
    exit()

print("4 - Lecture du fichier Fasta pour les ORF")
#lecture du fasta
if os.path.exists(fichier_orf):
    for record in SeqIO.parse(fichier_orf, "fasta"):
        desc = record.description
        match = re.search(r'location=(complement\()?(\d+)\.\.(\d+)', desc)  # Récupère les positions ORF
        if match:
            orf_start = int(match.group(2))
            orf_end = int(match.group(3))
            orf_name = re.search(r'\[protein=(.+?)\]', desc)  # Récupère le nom de l'ORF
            orf_name = orf_name.group(1) if orf_name else "Unknown_ORF"
            orf_data.append({"name": orf_name, "start": orf_start, "end": orf_end})
        else:
            print(f"Erreur : Position non trouvée pour {desc}")
    print("Lecture des ORF terminée !")
else:
    print("Fichier ORF.fasta introuvable. Aucune vérification des ORF ne sera effectuée.")


##Partie 5 : chercher si les variants sont sur des ORF

print("5 - Recherche d'ORF dans les variants")
def check_orf(mutations, orf_data):
    for mutation in mutations:
        mutation["ORF"] = "None"  #par défaut, aucune association
        #limites de la mutation
        mutation_start = mutation["pos"]
        mutation_end = mutation["pos"] + abs(mutation["svlen"])
        for orf in orf_data:
            orf_start = orf["start"]
            orf_end = orf["end"]
            #chevauchement
            if mutation_end >= orf_start and mutation_start <= orf_end:
                mutation["ORF"] = orf["name"]
                break  #dès qu'un ORF est trouvé on arrête la boucle
    return mutations

mutations_apparues = check_orf(mutations_apparues, orf_data)
mutations_disparues = check_orf(mutations_disparues, orf_data)
mutations_conservees = check_orf(mutations_conservees, orf_data)

##Partie 6 : sauvegarde des résultats
print("6 - Mise à jour des fichiers de résultats avec les ORF")

#enregistrement en fonction de la condition
comparaison_dir = os.path.join(dossier_filtrage, "comparaison_ORF_mutation")
if not os.path.exists(comparaison_dir):
    os.makedirs(comparaison_dir)
print(f"Les résultats seront enregistrés dans {comparaison_dir}")    

df_apparues_cold = pd.DataFrame([m for m in mutations_apparues if m["condition"] == "cold"])
df_apparues_hot = pd.DataFrame([m for m in mutations_apparues if m["condition"] == "hot"])
df_disparues_cold = pd.DataFrame([m for m in mutations_disparues if m["condition"] == "cold"])
df_disparues_hot = pd.DataFrame([m for m in mutations_disparues if m["condition"] == "hot"])
df_conservees_cold = pd.DataFrame([m for m in mutations_conservees if m["condition"] == "cold"])
df_conservees_hot = pd.DataFrame([m for m in mutations_conservees if m["condition"] == "hot"])

df_apparues_cold["Passage"] = [['P65']] * len(df_apparues_cold)
df_apparues_hot["Passage"] = [['P65']] * len(df_apparues_hot)
df_disparues_cold["Passage"] = [['P50']] * len(df_disparues_cold)
df_disparues_hot["Passage"] = [['P50']] * len(df_disparues_hot)
df_conservees_cold["Passage"] = [['P50', 'P65']] * len(df_conservees_cold)
df_conservees_hot["Passage"] = [['P50', 'P65']] * len(df_conservees_hot)

df_apparues_cold.to_csv(os.path.join(comparaison_dir, "mutations_apparues_cold_ORF.csv"), index=False, sep=",", quoting=1)
df_apparues_hot.to_csv(os.path.join(comparaison_dir, "mutations_apparues_hot_ORF.csv"), index=False, sep=",", quoting=1)
df_disparues_cold.to_csv(os.path.join(comparaison_dir, "mutations_disparues_cold_ORF.csv"), index=False, sep=",", quoting=1)
df_disparues_hot.to_csv(os.path.join(comparaison_dir, "mutations_disparues_hot_ORF.csv"), index=False, sep=",", quoting=1)
df_conservees_cold.to_csv(os.path.join(comparaison_dir, "mutations_conservees_cold_ORF.csv"), index=False, sep=",", quoting=1)
df_conservees_hot.to_csv(os.path.join(comparaison_dir, "mutations_conservees_hot_ORF.csv"), index=False, sep=",", quoting=1)

print("Résultats mis à jour avec les ORF et enregistrés ! ")
