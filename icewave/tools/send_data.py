import subprocess
import phonefleet.server.data as dataphone
import numpy as np
import phonefleet.network.whoami as whoami
import os
from datetime import datetime
from pathlib import Path

global LOCAL_FOLDER
LOCAL_FOLDER = '/data/data/com.termux/files/home/storage/downloads/Gobannos/'
global DATABASE
DATABASE = '/data/data/com.termux/files/home/storage/downloads/Gobannos/Database_serveur.txt'
#from LAITA computer
global SERVEUR_FOLDER #Folder on the server to send the data to
SERVEUR_FOLDER = '/Users/stephane/CASCADE26/Data/Buffer/Tests/'
global USERNAME #Username on the server to send the data to
USERNAME = 'stephane'
global IP_ADDRESS #IP Adress of the server to send the data to
IP_ADDRESS = '172.30.0.5'

#Create a database on phone to keep track of sended files. 
#Only send latest data that has been unmodified for more than one minute
import argparse


def gen_parser():   
    parser = argparse.ArgumentParser()
    parser.add_argument("-full", default = False,action="store_true",dest="full", help="Send all previous complete data set")
    args = parser.parse_args()
    return args

def command(filename):
    whoami_info = whoami.get_phone_info()
    print(f"Phone info: {whoami_info}")
    folder = f"{SERVEUR_FOLDER}{whoami_info['Name']}/" #also add date to the folder path
    filedest = '.'.join(os.path.basename(filename).split('.')[:-1])+ f'_{whoami_info['Name']}.csv' #remove the extension and name of the phone in the name of the files
    destination = f"{folder}{filedest}"
    cmd  = [
        "scp",
        "-i", "/data/data/com.termux/files/home/.ssh/id_ed25519",  # Assurez-vous que ce chemin est correct
        filename,
        f"{USERNAME}@{IP_ADDRESS}:{destination}"
    ]
    cmd = [
        "rsync",
        "-av",
        "--ignore-existing",
        "-e", "ssh -i /data/data/com.termux/files/home/.ssh/id_ed25519",
        f"--rsync-path=mkdir -p {folder} && rsync", # Creates all folders first
        filename,
        f"{USERNAME}@{IP_ADDRESS}:{destination}"
    ]
    #if not os.path.exists(folder):
    #    print(f"Create folder {folder} on the server !")#    os.makedirs(folder)
    return cmd

def send_data(data):
    for filename in data['filelist']:
        if has_been_sended(filename):
            print(f"{filename} has already been sent, skipping")
            continue
        print(f"Sending {filename} to the server...")
        try:
            #run a scp per filename, using the command function to build the command    
            result = subprocess.run(command(filename), check=True, capture_output=True, text=True)
            print("Fichiers copiés avec succès !")
            data_sended(filename)
        except subprocess.CalledProcessError as e:
            print(f"Erreur scp (Code {e.returncode}) : {e.stderr}")

def data_sended(filename):        
    with open(DATABASE, 'a') as f:
        f.write(f"{filename}\tsended\n")

def has_been_sended(filename):
    if os.path.exists(DATABASE):
        with open(DATABASE, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(filename) and 'sended' in line:
                    return True
    return False

def elapsed_time_since_modification(filename):
    chemin_fichier = Path(filename)
    if not chemin_fichier.exists():
        print(f"Le fichier {filename} n'existe pas.")
        return None

    # Étape 1 : Récupérer la date de modification
    timestamp = chemin_fichier.stat().st_mtime
    date_modification = datetime.fromtimestamp(timestamp)

    # Étape 2 : Calculer la différence avec l'heure actuelle
    temps_ecoule = datetime.now() - date_modification
    return temps_ecoule

def last_complete_file(data):
    sorted = dataphone.sort_data(data)
    complete = dataphone.complete_data(sorted)
    if len(complete)>0:
        time_elapsed = {}
        for num in complete:#only send complete data set
            elapsed = []
            for filename in sorted[num]['filelist']:
                elapsed.append(elapsed_time_since_modification(filename))
                #print(f"File {filename} last modified {elapsed[-1]} ago")
            time_elapsed[num] = np.mean(elapsed)

        if np.min(list(time_elapsed.values()))>np.timedelta64(60,'s'):
            num = complete[np.argmin(list(time_elapsed.values()))]
        else:
            print("Latest complete set of files has been modified recently, please wait")
            return None

        #useless to load the data !
        #sorted[num].update(dataphone.combine(sorted[num]))
        return sorted[num]
    else:
        print("No complete unmodified data found")
        return None


if __name__ == "__main__":
    args = gen_parser()
    
    datas = dataphone.search()
    keys = np.sort(list(datas.keys()))

    if args.full == True:
        for key in keys:
            data = datas[key]
            sorted = dataphone.sort_data(data)
            complete = dataphone.complete_data(sorted)
            for num in sorted.keys():
                if num in complete:
                    print(f"Sending complete data set {num} from {key}")
                    data_to_send = sorted[num]
                    send_data(data_to_send)
    else:
        last = keys[-1]
        data = datas[last]
        #sorted = dataphone.sort_data(data)
        #complete = dataphone.complete_data(sorted)
        data_to_send = last_complete_file(data)
        send_data(data_to_send)