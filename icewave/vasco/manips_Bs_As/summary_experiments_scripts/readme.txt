étape d'execution des codes :


lancer add_epsilonc_epsilondot_in_forcedisp_dicts.py si nouvelles données ajoutées 
pour les fit force vs disp 
et add_temperatures_in_sigmac_dicts_ice si nécessaire (si nouvelles données sigmac)

d'abord, séparément, exécuter 

group_all_force_displacement_data.py
et 
group_all_sigmac_data.py

(si des nouvelles données ont été ajoutées : on a donc exécuté au préalable 
le code "force_vs_displacement_from_paramsallacq_2cams_v2.py" qui est dans 
"D:/manips_BsAs/Summary/tracking_force_displacement/params_acquisitions/")

Ensuite il faudra aussi exécuter le code add_epsilonc_epsilondot_in_forcedisp_dicts.py
(pour ajouter les mesures d'epsilonc et epsilondot, pas faites avant ça)

Si jamais informations sur la temperature manquent (normalement c'st bon), on peut
lancer aussi add_temperatures_in_sigmac_dicts_ice.py

Ensuite on peut grouper les données qui matchent pour sigmac et mesures de
 modules elastiques effectifs avec :
group_sigmac_yougmodulus_data.py

Puis, enfin le code qui sert à lire et afficher les résultats est open_all_sigmac_force_disp_matching_data.py