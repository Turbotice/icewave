#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle

import csv
from datetime import datetime
import pytz
import glob 
import scipy
import sys
%matplotlib qt
# %%
csv_file_riki = 'D:/copie_BicWin25_geophones/a_copier_sur_backup25/General/Summary_geophone_lines/water_level_rimouski/02985_data.csv'

with open(csv_file_riki, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data = []
    for lines in csvFile:
        if count==6:
            header = lines
        elif count<6:
            pass
        elif count>=7:
            data.append(lines)
        count+=1

data = np.array(data)

waterlevels=np.zeros(np.shape(data)[0])

for i in range(len(data[:,0])):
    print(data[i,1])
    waterlevels[i] = float(data[i,1])

print(data)    

# %%
plt.figure()
plt.plot(waterlevels)
plt.show()
# %% conversion dates au format souhaité

def convert_date(date_str: str, tz_offset: str = "-0500") -> str:
    """
    Convertit une date du format 'YYYY/MM/DD HH:MM' 
    vers le format 'YYYY-MM-DDTHH:MM:SS.ssssss±HHMM'.

    Args:
        date_str (str): date au format 'YYYY/MM/DD HH:MM'
        tz_offset (str): décalage horaire au format ±HHMM (par défaut -0500)

    Returns:
        str: date convertie
    """
    # On parse la date d'entrée
    dt = datetime.strptime(date_str, "%Y/%m/%d %H:%M")
    
    # On reconstruit au format demandé
    formatted = dt.strftime(f"%Y-%m-%dT%H:%M:%S.%f{tz_offset}")
    
    return formatted


# Exemple d'utilisation
print(convert_date("2025/04/30 03:01"))

dates_newformat = []

for i in range(len(data[:,0])):
    #print(data[i,0])
    dtnfmt = convert_date(data[i,0])
    #print(dtnfmt)
    dates_newformat.append(dtnfmt)

# %%
np.mean(waterlevels)
#H_riki = 2.288541075831414 # approx bathymetrie
tidelevels = waterlevels - np.mean(waterlevels)
plt.figure()
plt.plot(tidelevels)
plt.show()
# %% def savetides function
import os
import csv
from datetime import datetime

def save_tides(dates, tidelevels, base_path, mmdd, overwrite=False):
    """
    Sauvegarde les données de marée pour un jour donné dans un fichier CSV structuré.

    Args:
        dates (list[str]): dates au format 'YYYY-MM-DDTHH:MM:SS.ssssss-0500'
        tidelevels (list[float]): niveaux de marée correspondants
        base_path (str): chemin de base, ex: 'D:/copie_BicWin25_geophones/Data'
        mmdd (str): date sous la forme 'mmdd' (ex: '0218')
        overwrite (bool): si True, écrase le fichier CSV existant

    Raises:
        FileExistsError: si le fichier existe déjà et overwrite=False
        ValueError: si aucune donnée ne correspond à la date choisie
    """

    # Filtrer les données correspondant au mmdd choisi
    selected_dates = []
    selected_tides = []
    for d, t in zip(dates, tidelevels):
        dt = datetime.strptime(d[:10], "%Y-%m-%d")  # extraire YYYY-MM-DD
        if dt.strftime("%m%d") == mmdd:
            selected_dates.append(d)
            selected_tides.append(t)

    if not selected_dates:
        raise ValueError(f"Aucune donnée trouvée pour la date {mmdd}")

    # Récupérer l'année et construire nom du fichier
    year = datetime.strptime(selected_dates[0][:10], "%Y-%m-%d").year
    filename = f"{year}-{mmdd[:2]}-{mmdd[2:]}_r01m_tides.csv"

    # Construire le chemin du dossier Marees
    target_dir = os.path.join(base_path, mmdd, "Marees")
    os.makedirs(target_dir, exist_ok=True)

    filepath = os.path.join(target_dir, filename)

    # Vérifier si le fichier existe déjà
    if os.path.exists(filepath) and not overwrite:
        raise FileExistsError(f"❌ Le fichier existe déjà : {filepath}")

    # Écriture du CSV (écrasement si overwrite=True)
    with open(filepath, mode="w", newline="") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(["date", "tide_height"])  # en-tête
        for d, t in zip(selected_dates, selected_tides):
            writer.writerow([d, f"{t:.3f}"])

    if overwrite:
        print(f"⚠️ Fichier écrasé et remplacé : {filepath}")
    else:
        print(f"✅ Fichier enregistré : {filepath}")


# %%
save_tides(dates=dates_newformat,tidelevels=tidelevels, base_path='D:/copie_BicWin25_geophones/Data/',mmdd='0304',overwrite=True)

# %%
