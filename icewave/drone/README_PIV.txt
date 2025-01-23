# Utilisation des Codes Drones/PIV_banquise

Ce fichier décrit l'utilisation des différents codes nécessaires au traitement d'une vidéo de Drone par le logiciel PIVlab sur Matlab.
Il ne décrit que les étapes permettant le traitement total des vidéos. 


## Passer de la vidéo '.MP4' à un dossier d'images '.tiff' : 
--> convert_multivideo2tiff.py : 
	Code Python, converti une série de vidéos en un seul dossier contenant l'ensemble des images au format : 'im_000.tiff'

## Traiter les images avec PIVlab : 

  ### Avant-propos : 
  --> La Correlation d'Images Digitales (DIC) nécessite l'entrée de plusieurs paramètres : 
	Dt : pas de temps entre deux images permettant de calculer le champ de vitesse
	W : taille de fenêtre au sein de laquelle les pixels sont comparés entre l'image i et l'image i+Dt
	b : pas de temps entre deux comparaisons d'images : u[i] = img[i+1] - img[i] , u[i+1] = img[i+b+1] - img[i+b] 

  --> Règles permettant une bonne PIV : 
  Sur l'ensemble de la zone d'intérêt, le champ de vitesse calculé doit satisfaire les conditions suivantes : 
	- u > 0.1 pix/frame => distinction du signal et du bruit
	- u < W/4 => le déplacement des pixels ne doit pas être trop grand par rapport à la taille de chaque boîte
	- u < ?? Eviter une trop grande rotation des champs de vitesse...

  ### Déterminer le pas de temps Dt nécessaire entre deux images : 
  --> Ouvrir PIVlab, importer des images séparées de [0,1,2,3,4,5,6,7...] images. Exécuter l'analyse PIVlab avec les paramètres souhaités,
	pre-processing des images, nombre de passage etc... 
  --> Balayer les champs de vitesse obtenus, sélectionner le champ de vitesse permettant d'avoir un signal suffisamment fort 
	dans la zone d'intérêt : u > 0.1 pix/frame. 
	Cela permet de distinguer le signal du bruit environnant. 
	Le champ de vitesse ne doit pas être trop élevé u < W/4, sinon la taille de fenêtre est trop faible pour le pas de temps Dt choisi... 
  --> Idéal : choisir Dt le plus faible possible, permettant d'assurer u > 0.1 pix/frame pour l'ensemble de la zone étudiée. 

  --> Cette étape peut être répétée sur plusieurs séquences d'images !! 

  ### Réaliser la PIV
  --> PIV_processing/automaticcycle_banquise_drone.m
	Permet de réaliser la PIV en utilisant PIVlab en parallel-computing
	Bien prendre garde à sélectionner le fichier que l'on souhaite traiter. Donc exécuter la première cellule du code avant de passer à la suite. 
	Le code créer un fichier .mat contenant : le champ de vitesse (u,v), la position des boîtes de PIV en pixel ainsi que les tables de paramètres (s et p) utilisées 
	pour effectuer la PIV .
	Le code appelle la fonction PIVlab_commandline_parallel_stereo décrite ci-dessous


  --> PIVlab_commandline_parallel_stereo.m
	Cette fonction prend plusieurs arguments :
	 - directory = fichier où se trouve les images .tiff à traiter (le nom complet du fichier doit être renseigné)
	 - reference = image de référence si l'on souhaite effectuer la PIV à partir d'une même image de référence, si égal à '', alors on compare les images 
		séparées de Dt deux à deux
	 - N = nombre de frames à process, si N = 0 alors la fonction process l'ensemble des images se trouvant dans le directory 
	 - i0 = première image à partir de laquelle on réalise la PIV
	 - Dt = pas de temps décrit précédemment
	 - b = décrit précédemment
	 - s = table de paramètres utiles à la réalisation de la PIV
	 - p = table de paramètres utiles au pre-processing des images 

	 La fonction appelle piv_analysis.m 
	 
	 Puis on réalise une boucle sur l'ensemble des images, les images sont preprocess à l'aide de la fonction PIVlab_preproc qui prend en argument l'image et la table de paramètres p. 
	 - PIVlab_preproc est une fonction implémentée dans PIVlab, permettant de préprocess une image à partir d'une table de paramètres donnée (cf documentation de PIVlab_preproc). 
	 Dans notre cas, nous utilisons uniquement une accentuation de contraste en utilisant un algorithme de CLAHE (Contrast Enhancement). 
	 Le contraste est alors localement accentué au sein d'une fenêtre de taille définie dans la table de paramètres p 'CLAHE size'. 
	 
	 La fonction appelle piv_FFTmulti (décrite ci-dessous) grâce à laquelle on récupère le champ de vitesse (u0,v0) entre les deux images comparées. 
	 - piv_FFTmulti est une fonction implémentée dans PIVlab, elle permet de réaliser la PIV et d'obtenir le champ de vitesse associé après corrélation entre deux images.
	 
	 
  --> piv_analysis.m 
	Permet de réaliser la PIV sur une paire d'images (filename1, filename2). Les actions réalisées par la fonction sont entièrement décrites dans la fonction Matlab. 
	Cette fonction peut être utilisée à l'aide du script 'PIV_image_pair.m'
	

  ### Post-processing des champs de vitesse brute 
  --> Main_data_structuration.m
	Ce script permet de traiter et ordonner les données brutes obtenues à partir de PIVlab. Les données brutes sont traitées par la fonction 'PIV_banquise_postprocessing.m' 
	décrite ci-dessous.
	Les données traitées ainsi que les paramètres de PIV sont ensuite assemblés au sein d'une structure matlab, grâce à la fonction 'genere_structure_banquise.m'. 
	Différents champs associés aux paramètres du drone, horaires et positions GPS sont ajoutés à cette structure, formant la structure '*_total_processed.mat'

	A partir de cette première structure, il est possible de générer une nouvelle structure à l'aide des fonctions 'scaling_structure.m' et 'scaling_structure_oblique.m'
	Ces fonctions permettent de mettre à l'échelle les champs de déplacement mesurés (y compris les champs de déplacement verticaux dans l'hypothèse où la surface filmée 
	est plane). A l'issue de l'une de ces fonctions, une structure matlab '*_scaled.mat' est sauvegardée. 

	Enfin l'ensemble des paramètres utilisés pour traiter les données est sauvegardé dans un fichier .txt à l'aide de la fonction 'saving_parameters.m'


  --> PIV_banquise_postprocessing.m
	Réalise le post-process des données brutes, obtenues à partir de PIVlab_commandline_parallel_stereo.m. Ce script prend comme argument :
	- le champ de vitesse brut (u,v) calculé par la PIV et les tables de paramètres utilisées (s,p)
	- W, la fenêtre utilisée par la PIV en pixels (dernière taille de la zone d'interrogation utilisée)
	- N, le nombre de frames à traiter, si N = 0, alors il faut traiter toutes le frames

	Il renvoie un champ de vitesse (u_filt,v_filt) filtré. 

	Le script applique une filtre médian à l'ensemble du champ de vitesse, retire les données abérrantes puis réalise une interpolation afin de remplacer les données abbérantes par des valeurs plus plausibles
	La méthode de filtrage utilisée par ce script est décrite dans les articles scientifiques suivants : 'Universal outlier detection for PIV data' Westerweel & Scarano (2005)


  --> genere_structure_banquise.m 

	Fonction permettant de créer une structure Matlab (un dictionnaire) à partir d'un fichier .mat généré par le code automaticcycle_banquise_drone.m. 
	Les arguments de cette fonction sont convenablement définis dans la fonction. 
	
	Elle retourne une structure Matlab contenant 
	- le champ de vitesse (m.Vx , m.Vy) dimensions [nx,ny,nt] (non scalés)
	- la grille spatiale utilisée pour la PIV (m.x , m.y) 
	- d'autres paramètres...
	
	Il est possible de supprimer les boîtes de PIV situées sur le contour de l'image étudiée en jouant sur le paramètre a. Si a = 1, alors toutes les boîts 
	situées sur le contour sont supprimées. Si a = 2, le contour formé par les deux séries de boîtes de PIV les plus à l'extérieures sont supprimées. 
	Cela peut être utile pour avoir toujours la même taille de boîte de PIV (en pixel). 
	
	
## Analyse des données après DIC/PIV

  ### Avant-propos 
  --> Les champs de vitesse obtenus après la PIV peuvent être analysés de différentes façons, les fonctions permettant d'obtenir les principales caractéristiques 
  du champ de vitesse sont décrites ci-dessous

  
  ### Main functions
  --> get_histogram_displacement.m :
  Fonction permettant d'obtenir l'histogram des déplacements moyen (moyenné en temps) pour chaque boîte de PIV. Donc de vérifier les critères de validation de la PIV

  --> plot_located_profile.m : 
  Permet d'obtenir le signal temporel d'une boîte PIV au cours du temps

  --> supress_quadratic_noise.m : 
  Supprime les variations quadratiques du champ de vitesse, associées aux variations de position du drone

  --> plot_velocity_features.m : 
  Plot le champ moyen, et écart-type du champ de vitesse pour Vx et Vy

  --> movie_velocity_field.m :
  Crée un film du champ de vitesse, en fixant la colorbar avec caxis_amp[0] et caxis_amp[1], affiche le film  avec un frame rate contrôllé par la variable 'fps'

  --> temporal_FFT.m : 
  Calcul la transformée de Fourier temporelle du champ de vitesse 

  --> plot_demodulated_field.m : 
  Plot le champ démodulé à différentes fréquences. Possible de créer un film et de sauvegarder les champs démodulés. 

  --> get_wave_vectors.m : 
  Calcul, pour chaque fréquence f, l'amplitude du vecteur d'onde k associé au champ démodulé à une fréquence f

  --> space_time_spectrum_Efk.m : 
  Calcul le spectre spatio-temporel E(f,k) (moyenné radialement dans l'espace des nombres d'onde) d'un champ de vitesse V

  --> extract_branchs_from_Efk.m : 
  Extraits les coordonnées (f,k) associées aux pics de plus hautes amplitudes dans le spectre spatio-temporel E(f,k)

  --> attenuation_coeff_corrected_direction.m :
  Calcul le coeffcient d'atténuation spatial pour différentes fréquences temporelles composantes d'un champ de vitesse V




