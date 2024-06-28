Extraction des datas issus des géophones (ISTerre) :

code ordinateur : pmmh

###################################
## Si problème / oubli => dossier SmartSolo, vidéo explicatives disponibles 
###################################

- Exécuter le logiciel Sololite.exe en tant qu'administrateur

#### Téléverser les données : 

- Placer les géophones dans l'ordre 1 -> 2 -> 3

- Créer un dossier correspondant à la date étudiée, avant d'exécuter le programme. Par exemple 0301/Geophones/

- Exécuter le script Python : move_data.py, qui permet de téléverser les fichiers DigiSolo vers l'ordinateur 

#########
Il faut bien relancer le code move_data.py tous les 3 géophones !!
#########

- Téléverser les fichiers log file dans /Canada_2024/

- Disconnect IGU (sur l'application Sololite)

- Une fois les lumières des géophones éteintes, on peut retirer les géophones et téléverser les géophones suivants (4 -> 5 -> 6) etc...

- Une fois tous les géophones downloadé : 'Export Seismic Data', rond bleu sur Sololite

- Choisir 'all components separate file'

- Changer les paramètres suivants : output data: miniSEED, Sample interval : 1 ms, s'assurer que cela correspond bien à facq, 

- Changer start_time et end_time associés aux fichiers du jour. Tous les fichiers du jour doivent se trouver dans la fenêtre de temps choisie

- Effectuer la séquence suivante sur Sololite : -> prepare, -> Vérifier le nombre de voies downloadées, -> Run ! 

- Les données téléversées se trouvent alors dans 'C:/SOLODATA/BicWin2024_12dB/BicWin2024_12dB/

- Utiliser le code prepare_solo.py afin de transformer les données, en utilisant la table de correspondance

#### Formattage des données 

- Le formatage des géophones se fait toujours 3 par 3 (1 -> 2 -> 3)

- Une fois les 3 géophones déposés, cliquer sur Setup IGU-16 = icône "carte SD"

- Effectuer les commandes suivantes : select all -> format device: IGU 16HR-3C -> write script -> Apply 

- Attendre que la lumière des géophones ne clignotte plus, la lumière rouge doit s'éteindre

- Ejecter les géophones depuis l'application, puis répéter l'opération avec la séquence de géophones suivante (3 -> 4 -> 5), etc... 