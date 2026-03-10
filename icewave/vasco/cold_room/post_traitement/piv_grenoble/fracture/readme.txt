le code le plus à jour faisant la reconstruction oblique et qui traite l'exemple du 29/11/2024 est piv_fracture_reconstruction_oblique_2.py

le code piv_fracture.py est plus basique, il prend simplement des echelles mesurées à la règle à 3 distances de la caméra et interpole pour avoir les rapports d'aspect
Il traite aussi l'exemple du 29/11.

si les echelles spatiales sont mesurées à plusieurs endroits sur la surface de la glace, le code piv_fracture_load_echelles_fracture.py charge le fichier où sont enregistrées ces 
échelles et interpole en 2D. Ensuite ça traite aussi l'exemple du 29/11