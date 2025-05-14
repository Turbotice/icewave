import os
import pandas as pd
import csv
import matplotlib.pyplot as plt 
import numpy as np
from victor_func import *
import threading
import pickle
from matplotlib.backend_bases import KeyEvent

phone_spacing = 3

nom_mesure = input("Nom de la mesure: ")


acqu_numb = input('acqu_numb: ')
if acqu_numb == '':
    acqu_numb = '1'
channel = input("Choisir la voie: \na_N \na_E \na_Z\n: ")
if channel == '':
    channel = 'a_Z'

base  = '../Data'
# Dictionnaire final
filename = '/pic_time/t0_time_' + nom_mesure  + '.pkl'
file2save = base + filename
if os.path.isfile(file2save):
    print('Time dictionnary already exists')
    with open(file2save,'rb') as f:
        time_dict = pickle.load(f)
else: 
    time_dict = {'d' + nom_mesure + 'a' + acqu_numb + 'tS' + '101' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '102' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '103' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '104' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '105' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '106' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '107' + 'N': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '101' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '102' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '103' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '104' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '105' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '106' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '107' + 'E': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '101' + 'Z': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '102' + 'Z': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '103' + 'Z': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '104' + 'Z': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '105' + 'Z': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '106' + 'Z': None,\
            'd' + nom_mesure + 'a' + acqu_numb + 'tS' + '107' + 'Z': None}
data = {}
dossier_csv = base +"/" + nom_mesure
# Parcours des fichiers du dossier
keys_channel = ['ts','a_N', 'a_E', 'a_Z',  'colonne5']
print("Chargement des données")

for fichier in os.listdir(dossier_csv):

    if fichier.endswith(".csv"):  # Vérifie que c'est un fichier CSV
        chemin_fichier = os.path.join(dossier_csv, fichier)

        cle = fichier.split('-')[0]
        # Chargement du fichier CSV en DataFrame
        df = pd.read_csv(chemin_fichier)
        # Stockage des colonnes dans un dictionnaire
        data[cle] = {keys_channel[i]: df.iloc[:, i].tolist() for i in range(df.shape[1])}


print("✅ Chargement terminé !")



directory = base + "/Tsync/"  
list_files(directory)
file_name= input("Choisir un fichier de syncro: ")

file_path = directory  + file_name
data_sync_tab = extract_data(file_path)
data_sync = {}
for ligne in data_sync_tab:

    cle = ligne[0].split('_')[0]
    if len(cle) == 1:
        cle = 'P0' + cle
    else:
        cle = 'P' + cle
    valeurs = ligne[1]  # Le reste des colonnes comme valeurs (liste)
    
    # Ajouter au dictionnaire
    data_sync[cle] = valeurs





data = dict(sorted(data.items()))
for i,phone_key in enumerate(data.keys()):
    data[phone_key][channel] = np.array(data[phone_key][channel], dtype = float)
    data[phone_key]['ts'] = np.array(data[phone_key]['ts'], dtype = float)
    data[phone_key]['ts'] = data[phone_key]['ts'][~np.isnan(data[phone_key][channel])]
    data[phone_key][channel] = data[phone_key][channel][~np.isnan(data[phone_key][channel])]


for cle in data.keys():
    data[cle]['ts'] = np.array(data[cle]['ts'], dtype = float)*10**-6  + np.array(data_sync[cle], dtype = float)
    print(data[cle]['ts'])

liste_tel = None  
def demander_input():
    global liste_tel
    liste_tel = selection_utilisateur(list(data.keys()))  
    plt.close()  
thread = threading.Thread(target=demander_input, daemon=True)
thread.start()

sample_rate = 1
fig, ax = plt.subplots()

for i,phone_key in enumerate(data.keys()):

    current_stream = data[phone_key][channel] - np.mean(data[phone_key][channel])
   

    ax.plot(data[phone_key]['ts'][::sample_rate],(current_stream/max(current_stream) + phone_spacing*(i+1))[::sample_rate], label = phone_key)
    
ax.set_title(f'Signal superposition for channel {channel}')
plt.legend()
plt.show()
thread.join()
is_zooming_or_panning = False

                   

selected_keys = list(time_dict.keys()) 
selected_keys = [element for element in  selected_keys if element[-1] == channel[-1]]
selected_index = 0 
vertical_lines = []
annotations = []
def update_annotations():
    """ Met à jour les annotations sans effacer le reste du graphique """
    global vertical_lines, annotations
    ax = plt.gca()

    # Supprimer les anciennes annotations
    for line in vertical_lines:
        line.remove()
    for text in annotations:
        text.remove()

    vertical_lines.clear()
    annotations.clear()

    # Ajouter les nouvelles annotations
    for key, value in time_dict.items():
        if value != None:
            line = ax.axvline(x=value, color='r', linestyle='--')
            text = ax.text(value, ax.get_ylim()[1] * 0.9, key[-5:-1], 
                       color='red', verticalalignment='bottom', fontsize=10)
            vertical_lines.append(line)
            annotations.append(text)

    fig.canvas.draw()  # Rafraîchir l'affichage

def on_click(event):
    global selected_index
    global selected_index, is_zooming_or_panning
    if is_zooming_or_panning:
        return  # Ignore la sélection si l'utilisateur est en train de zoomer ou déplacer

    if event.inaxes is not None and event.xdata is not None:  # Vérifie que le clic est dans la zone du graphique et qu'on n'a pas déjà 5 valeurs
        # Stocke la valeur dans le dictionnaire
        print("choisir une source:")
        liste_key_input = []
        for key in selected_keys:
            key_name = key[-5:-1]
            liste_key_input.append(key_name)
            print(key_name)
        key_input = input("Entrez la clé pour cette valeur sélectionnée: ")
        if key_input in liste_key_input:
            time_dict['d' + nom_mesure + 'a' + acqu_numb + 't' + key_input + channel[-1]] = event.xdata
            update_annotations()
        else:
            return



def on_pan_or_zoom(event):
    global is_zooming_or_panning
    is_zooming_or_panning = True 

def on_release(event):
    global is_zooming_or_panning
    is_zooming_or_panning = False


fig, ax = plt.subplots()

for i,phone_key in enumerate(liste_tel):

    current_stream = data[phone_key][channel] - np.mean(data[phone_key][channel])

    ax.plot(data[phone_key]['ts'][::sample_rate],(current_stream/max(current_stream) + phone_spacing*(i+1))[::sample_rate], label = phone_key)
    
ax.set_title(f'Signal superposition for channel {channel}\n Selectioner les sources')
fig.canvas.mpl_connect("button_press_event", on_click)
fig.canvas.mpl_connect("motion_notify_event", on_pan_or_zoom)  
fig.canvas.mpl_connect("button_release_event", on_release)  
plt.legend()
plt.show()

print("\n📌 Valeurs sélectionnées :")

time_dict_without_none = {k: v for k, v in time_dict.items() if v is not None}
for key, value in time_dict_without_none.items():
    print(key[-5:], ":", value)



with open(file2save, 'wb') as f:
    pickle.dump(time_dict, f)