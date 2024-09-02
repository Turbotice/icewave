import glob
import os
import numpy as np
import csv

global table
table = {'Accelerometer':'a','Gyroscope':'g','Location':'l','Magnetometer':'m'}
coords = {'x','y','z'}

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
       #try:
        print(key)
        print(datafile)
        data[key]={}
        data[key]['d'] = np.loadtxt(datafile, delimiter=',',usecols=(1,2,3),skiprows=1)
            #get times
        data[key]['t'] = np.loadtxt(datafile, delimiter=',',dtype=str,usecols=(0),skiprows=1)
        #except:
            #print('data cannotbe converted to float')
    return data

def sort(data):
    newdata = {}
    keys = {'time','device'}
    for key in keys:
        if key in data.keys():
            newdata[key]=data[key]
    
    keys = {'a','m','g'}
    
    for key in keys:
        if key in data.keys():
            for i,k in enumerate(coords):
                newdata[key+k]=data[key]['d'][:,i]
        newdata['t'+key] = data[key]['t'].astype(float)

    key = 'l'
    if key in data.keys():
        newdata['loc'] = {}
        newdata['loc']['lat']=data[key]['d'][:,0]
        newdata['loc']['lon']=data[key]['d'][:,1]
        newdata['loc']['elev']=data[key]['d'][:,2]
        newdata['loc']['t']= data[key]['t'].astype(float)
    
    data['coords']=coords
    return newdata

def read_meta(csvfile):
    d={}
    with open(csvfile) as f:
        csv_reader = csv.reader(f, delimiter=',')
        line=0
        for row in csv_reader:
            if line==0:
                headers = row
            else:
                #print(headers,row)
                for i,header in enumerate(headers):
                    d[header[:6]+'_'+row[0]]=row[i]
            line+=1

    return d

