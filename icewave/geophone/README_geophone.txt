Traitement de données Géophones ISTerre : 

Ce document README décrits les différentes étapes à suivre afin de déterminer le module d'Young, le coefficient de Poisson
ainsi que l'épaisseur de glace lors d'une mesure par une ligne de géophones


## Raw Datas : 

Après extraction des données par le logiciel associé aux géophones. On extracte les fichiers suivants : 
- un fichier miniseed pour chaque voie d'un géophone (Z,E,N) et cela pour chaque acquisition 
- un fichier .LOG associé à chaque géophone, qui nous donne le temps de début de chaque acquisition. Plusiseurs acquisitions
peuvent être enregistrées dans un même fichier .LOG

ATTENTION : les géophones sont nommés par leur numéro de série, à chaque numéro de série est associé un géophone. 
Pour renommer les fichiers .miniseed, il faut utiliser le code 'prepare_smatsolo.py' décrit ci-dessous. 


## Préparation des données géophones : 
--> prepare_smartsolo.py : code permettant de formater les données géophones. 

Dans ce code, il faut indiquer : 
- le path des données miniseed
- le path des fichiers .LOG associé
- le path de la table de correspondance

Ce code permet de trier les données géophones par acquisition et de renommer les fichiers .miniseed correctement

## Détermination des temps t0 :
--> extract_times.py : Permet d'afficher les signaux de tous les géophones sur chaque voie (E,N,Z) puis de déterminer 
le temps UTC auquel chaque source a été réalisée. On ne garde qu'une seule source (Z,E ou N) pour chaque point source. 
Chaque point source est défini par un matricule S101, S102, S103, S104, S105, S106

Ces temps sont sauvegardés dans un dictionnaire, chaque source étant définie par une clé du type:
'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 

## Relations de dispersion :
--> ZEN_dispersion.py : ce code permet d'obtenir les relations de dispersion associées à chaque mode (longitudinal QS0 (voie N),
cisaillement SH0 (voie E), flexion QS (voie Z)). Fonctionnement du code ci-dessous:

1) Load des signaux

Pour chaque type de source (E,N,Z) et chaque direction (1 ou 2)
2) Load des temps t0 de chaque point source associé à une direction et une composante (E,N,Z) = channel (0,1,2)
et construction d'une matrice des signaux
3) Calcul du spectre FK à l'aide d'une méthode SVD
4) Affichage du spectre FK et détermination des vitesses de phase pour les ondes longitudinales et de cisaillement (QS0, SH0)
ou extraction des coordonnées (f,k) pour le mode de flexion (QS)


## Détermination des propriétés mécaniques:
--> ice_properties_Enuh.py : inversion des propriétés mécanique de la glace (E,nu,h) à partir des vitesses de phase (c_QS0, c_SH0) et
de la relation de dispersion (f,k) du mode de flexion. Fonctionnement du script:

1) Load des vitesses de phase (c_QS0, c_SH0) et de la relation de dispersion du mode de flexion pour chaque direction (1 et 2)
2) Calcul du module d'Young et du coefficient de Poisson
3) Inversion de l'épaisseur de glace h en approchant la relation de dispersion (f,k)

