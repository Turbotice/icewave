
import glob
import os
#Global variables0
import zipfile as zip

import platform
import socket
import csv

import numpy as np
import icewave.phone.time_phone as time_phone

import icewave.tools.rw_data as rw

global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()


def find_path(base,disk='labshared2'):
    print('OS type : '+str(osname))

    if 'Windows' in ostype:    
        if disk=='labshared2':
            serveurfolder = 'W:/'+base #fucking windows OS : beware of the mounting disk you used
    
    if 'Linux' in ostype:
        if 'adour' in osname:
            serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base
        elif 'spi201711-Latitude-5480' in osname:
            serveurfolder = '/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/'+disk+'/'+base
        elif 'thiou' in osname:
            serveurfolder = '/volume3/'+disk+'/'+base #praise UNIX system    		
        else:
            serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base #praise UNIX system
			
    if 'Darwin' in ostype:# or 'laita' in osname:    
        serveurfolder = '/Volumes/'+disk+'/'+base #praise UNIX system        

    return serveurfolder	

def extract_all(folder):
#    folder = 'Telephones/Soufflerie_dec23/131223/Telephones/121223_4_U400cms/'
    filelist = get_ziplist(folder)
    print(filelist,folder)
    for zipfile in filelist:
        print(zipfile)
        extract(zipfile)

def get_phone_list(folder):
    flist = []
    folderlist = glob.glob(folder+'*/')
    for folder in folderlist:
        csvfiles = glob.glob(folder+'*.csv')
        if len(csvfiles)>0:
            #print(get_phone_num(folder))
            flist.append(folder)
    return flist
            
def get_ziplist(folder):
   # folder = find_path(path)
    return glob.glob(folder+'*.zip')

def extract(zipfile):
    with zip.ZipFile(zipfile, 'r') as zip_ref:
        extractfolder = zipfile.split('.zip')[0]
        zip_ref.extractall(extractfolder)

def get_phone_num(filename):
    return int(filename.split('_')[-1].split('/')[0])-100

def get_filename_key(filename):
    return os.path.basename(filename).split('.csv')[0]

def get_keys(name):
    print(name)
    if name == 'Accelerometer':
        return ['ta','ax','ay','az']
    if name == 'Magnetometer':
        return ['tm','mx','my','mz']
    if name == 'Gyroscope':
        return ['tg','gx','gy','gz']
    if name == 'Location':
        return ['tl','Lat','Long','H','V','Dir','Hacc','Vacc']
    if name == 'time':
    	return ['event','expt','syst','systext']
    if name == 'device':
        return ['property','value']

def load_data(folder):
    csvlist = glob.glob(folder+'*.csv')
    
    csvheader = glob.glob(folder+'*/*.csv')
    print(csvheader)

    data={}
    for csvfile in csvlist+csvheader:
        print(csvfile)
        d = load_csv(csvfile)
        data.update(d)
            
    return data

def load_pickle_data(folder):
    filelist = glob.glob(folder +'*/*phonedata.pkl')
    data = {}
    for filename in filelist:
        phone = filename.split('/')[-2].split('_')[1]
        print(phone,filename)
        d = rw.load_pkl(filename)
        data[phone]=d
    return data

def load_csv(filename):
    data = {}
    rows = []

    keys = get_keys(get_filename_key(filename))
    #print(keys)
    with open(filename) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for i,row in enumerate(spamreader):
            rows.append(row)
            if len(keys)>0:
		        #print(', '.join(row))
                for j,k in enumerate(keys):
                    if i==0:
                        data[k]=[]
                    else:
                        try:
                            s = float(row[j])
                        except:
                            s = row[j]
                        data[k].append(s)
            else:
                print('following variables not catched :')
                print(row)
    return data 

def convert_dict(data,phone):
    header = ['Phone']+list(data.keys())
    datamat = [phone]+[data[key] for key in data.keys()]
    return header,datamat

def convert_super_dict(results):
    #keys contain the phonelist
    datamat=[]
    phonelist = []

    phonelist = list(np.sort(list(results.keys())))
    for key in phonelist:
        header,data = convert_dict(results[key],key)
        #print(len(header))#check that the header have all the same length
        datamat.append(data)
    return header,datamat

def save_data_single_phone(data,savefolder,suffix=''):
    import pickle
    if 'time' in data.keys():
        found = True
        if 'system_START' in data['time'].keys():
            key = 'system_START'
        elif 'system time text_START' in data['time'].keys():
            key = 'system time text_START'
        else:
            found=False
            print('time START key not found')
            print('no date available')
            print(found)
        if found:    
            month = '0'+str(time_phone.get_time(data['time'][key])[1].month)
            day = str(time_phone.get_time(data['time'][key])[1].day)    
        else:
            month = ''
            day = ''
        date = month+day
    else:
        print('no date available')
    filename = savefolder + '_phonedata'+suffix+'.pkl'#+date+'_'+str(phone)+
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def write_csv(data,savefolder,title=''):
    filename = savefolder+title+'_Summary.csv'
    header,datamat = convert_super_dict(data)

    with open(filename, 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)

        # write multiple rows
        writer.writerows(datamat)
