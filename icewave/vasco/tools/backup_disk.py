#%%
import os
import shutil
#%%
def copier_dossier(source, destination):
    """
    Copie tout le contenu d'un dossier source vers un dossier destination,
    en évitant de recopier les fichiers ou dossiers déjà présents.
    """
    # Crée le dossier destination s'il n'existe pas
    if not os.path.exists(destination):
        os.makedirs(destination)

    # Parcours tous les fichiers et sous-dossiers dans source
    for root, dirs, files in os.walk(source):
        # Calculer le chemin relatif par rapport à la source
        relative_path = os.path.relpath(root, source)
        # Définir le chemin cible correspondant
        dest_path = os.path.join(destination, relative_path)

        # Crée le dossier dans la destination s'il n'existe pas
        if not os.path.exists(dest_path):
            os.makedirs(dest_path)

        # Copier les fichiers
        for file in files:
            source_file = os.path.join(root, file)
            dest_file = os.path.join(dest_path, file)
            # Copier seulement si le fichier n'existe pas déjà
            if not os.path.exists(dest_file):
                shutil.copy2(source_file, dest_file)  # copy2 conserve les métadonnées
                print(f"Copié : {source_file} -> {dest_file}")
            else:
                print(f"Déjà présent, ignoré : {dest_file}")
#%%
# Exemple d'utilisation
#source_folder = "D:/manips_BsAs"
#destination_folder = "E:/manips_BsAs"

# NE PAS SE TROMPER ET INTERVERTIR LES DEUX !!!!!

#copier_dossier(source_folder, destination_folder)