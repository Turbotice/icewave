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


## Affichage des signaux sur une même voie pour tous les géophones : 
--> plot_smartsolo.py : 

Ce code permet d'afficher une même voie pour chaque géophone afin de déterminer les signaux d'intérêts permettant de 
réaliser une inversion. 
Ce code contient plusieurs sections : 

-1-- Rentrée de variables : 
	path2data : path vers les données miniseed
	geophones_table_path : path vers la table de correspondance
	channel : voie des géophones que l'on veut afficher (Z, E ou N)

-2-- Identification des signaux : 
	Il faut d'abord identifier les signaux réalisé en chaque source (S1, S2, ...). Pour cela, on affiche les signaux 
sur la voie Z (voie sur laquelle on voit le mieux les signaux). 

Chaque bruit effectué devrait être mieux observable sur la voie associée : source N, plus visible sur la voie N, source E plus 
visible sur la voie E etc... Mais tous les signaux sont plus visibles sur la voie Z. 

Notations : sources S1, S2, S3 (proche de G1) direction dir1
	    sources S4, S5, S6 (proche de G16) direction dir2

	Il faut détecter, pour chaque type de signal (Z,E,N), trois signaux différents (un signal par source : S1, S2, S3, puis
S4, S5 et S6). On choisit un temps de départ t1 pour chacun de ces signaux. t1 correspond au temps auquel le signal créée arrive 
sur le premier géophone (G1 pour la direction dir1, G2 pour la direction dir2). 

On doit avoir, pour chaque voie (Z,E,N) un temps t1 par source, cela pour chaque direction : dir1 / dir2. 

-3-- FK-plot simple : 
	Une étape bonus est d'afficher un plot FK pour chaque signal envoyé. 
Cette étape nous permet de vérifier la qualité des signaux choisis. Si on observe correctement la relation de dispersion,
cela signifie que les signaux choisis sont de bonne qualité. 

-4-- Obtention des signaux pour une voie donnée : 
	Une fois les temps t1 déterminés pour une voie donnée, on peut sélectionner les signaux sur cette voie en modifiant
le channel en entrée du code 'plot_smart_solo.py'. Une fois le channel sélectionné, on peut vérifier que les signaux sélectionnés
sont bien visibles, et nous donne un FK simple. 

-5-- Création d'une matrice 3D : 
	On crée une matrice à partir des temps sélectionnés pour une seule voie donnée, pour une direction donnée. 
Cette matrice sera utilisée pour obtenir une relation de dispersion en utilisant une méthode en développement en valeur
discreète (VSD). 

-6-- Obtention du plot FK pour une voie donnée : 
	Pour obtenir ce plot FK pour les signaux d'une voie donnée pour tous les géophones, on utilise la matrice construite
précédemment : 'seismic_matrix', ainsi que le code python : 'FK_with_SVD.py'. 
Dans ce code, il faut préciser les distances entre géophones (en mètre), ainsi que le rang des valeurs singulières. 
On peut alors tracer la courbe FK, en utilisant la fonction fn_svd décrite dans le script. (voir article Ludovic Moreau)

-7-- Obtention des vitesses C_shear et C_longi : 
	A partir des plots FK pour les voies E et N, on peut retrouver les vitesses de propagation des ondes de cisaillement
(voie E), et des ondes de compression (voie N). Pour cela, on pointe sur chaque graphique FK la courbe affine associée à chaque 
type d'ondes. Cela se fait dans le code 'plot_smart_solo.py'

-8-- Obtention du module d'Young et du coefficient de Poisson : 
	A partir des vitesses C_shear et C_longi, on est capable de déterminer le module d'Young E, et le coefficient de Poisson
nu du matériau. code : 'plot_smart_solo.py'

-9-- Obtention de l'épaisseur de glace : 
	On peut enfin utiliser la voie Z et les ondes de flexion pour déterminer l'épaisseur de glace h. 
Pour cela, on récupère les signaux sur la voie Z (3 séries de signaux Z). On crée la matrice 3D, puis on utilise le code 
'FK_with_SVD.py' pour obtenir le plot FK des ondes de flexion. 
Sur ce plot, on peut alors pointer la relation de dispersion des ondes de flexion, ces points sont sauvegardés dans un fichier 
'dispersion_dir.pkl'. 
	En utilisant ces points, le module d'Young E et le coefficient de poisson nu, on est alors capable d'inverser la relation
de dispersion et d'obtenir l'épaisseur de glace en utilisant le code 'invert_FK.py'



 
	
