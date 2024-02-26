import glob
import os
import numpy as np
import csv

global table
table = {'Accelerometer':'a','Gyroscope':'g','Location':'l','Magnetometer':'m'}

def extract_all(folder):
    filelist = glob.glob(folder+'*.zip')
    print(f"Extract : {filelist}")

    for filename in filelist:
        with zipfile.ZipFile(filename,"r") as zip_ref:
            foldersave = filename.split('.')[0]
            zip_ref.extractall(foldersave)

def get_pĥonelist(folder):
    folders = glob.glob(folder+"*")
    phonelist = [get_number(f) for f in folders]
    return list(set(phonelist))

def get_number(folder):
    return int(os.path.basename(folder).split('_')[1])

def loads(folderlist):
    data = {}

    for folder in folderlist:
        phone = get_number(folder)
        print(f"Load data for phone {phone}")
        data[phone] = load(folder)
    return data

def load(folder):
    metalist = glob.glob(folder+'/meta/*.csv')
    datalist = glob.glob(folder+'/*.csv')

    data = {}
    for meta in metalist:
        key = os.path.basename(meta).split('.')[0]
        data[key]  = read_meta(meta)

    for datafile in datalist:
        key = os.path.basename(datafile).split('.')[0]
        if key in table:
            key = table[key]
        data[key] = np.loadtxt(datafile, delimiter=',',skiprows=1)
    return data

def read_meta(csvfile):
    d={}
    with open(csvfile) as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            d[row[0]]=row[1]
    return d

