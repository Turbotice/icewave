
import os
import csv
import numpy as np
def selection_utilisateur(liste):
    # Affichage des options
    print("\nSélectionnez des éléments en entrant un numéro ou un intervalle (ex: 1,3-5,7) : \n")
    print(liste)

    # Récupération de l'entrée utilisateur
    choix_utilisateur = input("\nVotre sélection : ")
    choix_valides = []

    # Traitement des entrées
    for morceau in choix_utilisateur.split(","):
        morceau = morceau.strip()  # Supprime les espaces inutiles

        if "-" in morceau:  # Gestion des intervalles (ex: "2-4")
            try:
                debut, fin = map(int, morceau.split("-"))
                if 1 <= debut <= fin <= len(liste):  # Vérifie que l'intervalle est valide
                    choix_valides.extend(range(debut - 1, fin))  # Ajoute l'intervalle
            except ValueError:
                print(f"⚠ Erreur : Intervalle invalide ({morceau})")

        else:  # Gestion des numéros simples (ex: "3")
            try:
                num = int(morceau)
                if 1 <= num <= len(liste):  # Vérifie que le numéro est valide
                    choix_valides.append(num - 1)
            except ValueError:
                print(f"⚠ Erreur : Numéro invalide ({morceau})")

    # Suppression des doublons et tri des choix
    choix_valides = sorted(set(choix_valides))

    # Récupération des éléments sélectionnés
    elements_choisis = [liste[i] for i in choix_valides]

    print(f"\n✅ Vous avez sélectionné : {elements_choisis}")
    return elements_choisis




def list_files(directory):

    for root, dirs, files in os.walk(directory):

        for file in sorted(files, key=lambda x: os.path.getmtime(os.path.join(root, x))):
            print(file)

def extract_data(file_path):
    data = []

    with open(file_path, "r") as file:
        reader = csv.reader(file)
        next(reader)  # Ignorer la première ligne (en-tête)
        
        for row in reader:
            if len(row) >= 2:
                num = row[0]  # Colonne '#'
                tlag = row[3] # Colonne 'tlag'
                data.append((num, tlag))
    
    return data

def detect_impulse(signal, sample_rate, window_size=2, step_size=0.05):
    window_samples = int(window_size * sample_rate)
    step_samples = int(step_size * sample_rate)
    max_energy = 0
    best_start = 0

    for start in range(0, len(signal) - window_samples, step_samples):
        window = signal[start:start + window_samples]
        energy = np.sum(np.square(window))

        if energy > max_energy:
            max_energy = energy
            best_start = start

    return best_start, best_start + window_samples

def choisir_sources(sources_disponibles):
    print("Sources disponibles :")
    for i, source in enumerate(sources_disponibles, 1):
        print(f"{i}. {source}")

    choix = input("Entrez les numéros des sources à utiliser (séparés par des virgules) : ")
    indices = [int(i) - 1 for i in choix.split(",") if i.strip().isdigit()]
    
    return [sources_disponibles[i] for i in indices if 0 <= i < len(sources_disponibles)]
