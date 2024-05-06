# Détails fichiers issus du Drone DJI Mavic 3 Pro

## Les différents types de fichier

- Vidéo au format .MP4
- Vidéo au format .LRF (Low Resolution File), peut être converti en .MP4 pour avoir un film de faible résolution 
- fichier format .SRT (SubRip Subtitle) file, utilisé pour ajouter des sous-titres aux vidéos
- photo au format .JPG 
- flightrecords => chercher comment ouvrir et récupérer ces documents ! 

## Organisation des fichiers : 

1. Trier les données par date, exemple : 2023/0130/ 
2. Indiquer le numéro d'expériences de la journée : 01,  02,  03 ...
3. Séparer les différents types de fichier entre chaque expérience : MP4, LRF, SRT, ...
4. Utiliser improc_serveur.py pour transformer les vidéos en dossier de fichier .TIFF
5. Faire de la PIV ou autre traitement 


## Récupération flightrecords
1. Connexion par câble USB avec la manette de commande 
2. Récupération des données de vol : " flightrecords" pour chaque vol effectué
3. Aller sur phantomHelp, load un fichier flightrecord 
4. Télécharger le fichier csv et le fichier KML associé
5. Placer les fichiers dans l'endroit associé