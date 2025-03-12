import matplotlib.pyplot as plt

# Fonction pour extraire les données de chaque section du VCF
def extract_vcf_data(vcf_filename):
    sections = {}
    with open(vcf_filename, 'r') as vcf:
        current_section = None
        for line in vcf:
            line = line.strip()  # Supprimer les espaces blancs en début et fin de ligne
            if line.startswith('##FILENAME'):
                # Nouveau fichier (section) trouvé
                current_section = line.split('=')[1].strip()
                sections[current_section] = []  # Créer une nouvelle entrée pour chaque section
            elif line.startswith('#CHROM') or not line:
                # Ignorer la ligne d'en-tête ou les lignes vides
                continue
            else:
                # Lignes de variantes
                parts = line.split('\t')
                # On attend au moins 10 colonnes (les 9 habituelles + TEMPERATURE)
                if len(parts) < 10:
                    continue  # Ignorer les lignes mal formatées
                chrom, pos, variant_id, ref, alt, qual, filter, info, format, *samples = parts
                # Extraire la température à partir du dernier champ
                temperature = samples[-1] if samples else None
                # Extraire la fréquence allélique (AF) du champ INFO
                allele_freq = None
                for field in info.split(';'):
                    if field.startswith('AF='):
                        allele_freq = float(field.split('=')[1])
                        break
                if allele_freq is not None:
                    # Stocker la position, l'allele frequency et la température
                    sections[current_section].append((pos, allele_freq, temperature))
    return sections

# Créer et enregistrer les graphiques
def plot_allele_frequencies(vcf_filename):
    sections = extract_vcf_data(vcf_filename)
    
    for section, variants in sections.items():
        if not variants:
            continue  # Si aucune variante n'a été trouvée pour cette section, passer à la suivante

        # Supposons que la température est identique pour tous les variants d'une section, on l'extrait du premier variant
        temp_val = variants[0][2] if variants[0][2] is not None else "No Temperature Info"

        # Créer une figure pour chaque section
        fig, ax = plt.subplots(figsize=(10, 5))
        
        positions = [v[0] for v in variants]
        allele_freqs = [v[1] for v in variants]
        
        ax.bar(positions, allele_freqs, color='skyblue')
        ax.set_title(f"Allele Frequencies for {section} - {temp_val}")
        ax.set_xlabel('Position')
        ax.set_ylabel('Allele Frequency')
        ax.tick_params(axis='x', rotation=45)
        
        # Sauvegarder chaque graphique en .png
        fig.tight_layout()
        plt.savefig(f"{section}_allele_frequencies.png")
        plt.close(fig)  # Fermer la figure pour libérer la mémoire

# Exemple d'utilisation avec votre fichier VCF
plot_allele_frequencies('combined.vcf')
