import numpy as np
import pylab as plt
import glob
import os
import shutil
from pprint import pprint

import fitdecode
try:
    import garmin
except:
    print('Please install garmin package')

import icewave.tools.datafolders as df

base = df.path_GPSdata()
print('')
print(base)

folderactivity = base+'GARMIN/Garmin/Activities/' # on macOS X
folderwpts = base+'GARMIN/Garmin/GPX/' # on macOS X

#destfolder = '/Users/stephane/Documents/Rimouski 2023/Data/GPS/'
#serveurfolder = '/Volumes/labshared2/Banquise/Rimouski 2023/Data/GPS/'

def ls_date(date=''):
    if date=='':
        date = df.get_current_date()
    dayGPS = date_GPS(date)

    fitlist = glob.glob(folderactivity+dayGPS+'*.fit')
    pprint(fitlist)
    
    wayptslist = glob.glob(folderwpts+'Way*'+dayGPS+'*')
    pprint(wayptslist)
    return fitlist,wayptslist
    
    
def date_GPS(date):
    year,day = date.split('_')
    print(year,day)
    month = day[:2]
    day = day[2:]
    return year+'-'+month+'-'+day

def download(date=''):
    fitlist,wayptslist = ls_date(date=date)

    base = df.instrument_folder(key='GPS',date=date)
    print(base)
    for i,filename in enumerate(fitlist):
        
        dest_file = base+str(df.ndigit(i,n=4))+'_'+os.path.basename(filename)
        dest_file = dest_file[:-4].replace('.','_')+'.fit'

        print(filename,dest_file)
        copy(filename,dest_file)

    for i,filename in enumerate(wayptslist):
        dest_file = base+str(df.ndigit(i+1000,n=4))+'_'+os.path.basename(filename)
        print(filename,dest_file)
        copy(filename,dest_file)


def download_all(date=None):
    """
    Load .fit and .gpx data from Garmin to destfolder & serveurfolder
    If date is not specified, look only for the data from the current date.
    INPUT 
    ----- 
    date : string 
    OUTPUT
    -----
    None
    """
    if date is None:
        date = df.get_current_date()
    else:
        filelist = glob.glob(folderactivity+'/'+date+'*')

    for filename in filelist:
        date = get_date(filename)
        print(date)

        dest_file1 = destfolder+date+'/'+os.path.basename(filename)
        print(dest_file1)
        copy(filename,dest_file1)
        
        dest_file2 = serveurfolder+date+'/'+os.path.basename(filename)
        print(dest_file2)
        copy(filename,dest_file2)
        

    filelist = glob.glob(folderwpts+'*.gpx')

    for filename in filelist:
        date = get_date_wpts(filename)
        print(date)
        
        dest_file1 = destfolder+date+'/'+os.path.basename(filename)
        print(dest_file1)
        copy(filename,dest_file1)

        dest_file2 = serveurfolder+date+'/'+os.path.basename(filename)
        print(dest_file2)
        copy(filename,dest_file2)
        
def copy(src_file,dest_file):
    destfolder = os.path.dirname(dest_file)
    if not os.path.exists(destfolder):
        os.makedirs(destfolder)
    shutil.copyfile(src_file, dest_file, follow_symlinks=True)

def get_date(filename):
    return filename.split('/')[-1].split(' ')[0].replace('-','')

def get_date_wpts(filename):
    return filename.split('/')[-1].split('_')[1].split('.')[0].replace('-','')
    
def decode(filename):
    d={}
    d['r']=[]
    with fitdecode.FitReader(filename) as fit:
        for frame in fit:
        # The yielded frame object is of one of the following types:
        # * fitdecode.FitHeader (FIT_FRAME_HEADER)
        # * fitdecode.FitDefinitionMessage (FIT_FRAME_DEFINITION)
        # * fitdecode.FitDataMessage (FIT_FRAME_DATA)
        # * fitdecode.FitCRC (FIT_FRAME_CRC)
            if frame.frame_type == fitdecode.FIT_FRAME_DATA:
            # Here, frame is a FitDataMessage object.
            # A FitDataMessage object contains decoded values that
            # are directly usable in your script logic.
                if frame.name=='record':
                    fsave = frame
                    d['r'].append(frame)
    return d

def convert_traj(d,tmin=None,tmax=None):
    """
    Convert trajectory from
    INPUT 
    ----- 
    d : dict with key 'r', with fields position_long and position_lat
    tmin : float
    tmax : float
    OUTPUT
    -----
    X,Y : np 1d arrays
    """
    X,Y = [],[]
    for i,a in enumerate(d['r']):
        x = a.get_value('position_long')
        y = a.get_value('position_lat')
        if x is not None:
            if tmin is not None:
                b1 = a.get_value('timestamp').time()>tmin
                b2 = a.get_value('timestamp').time()<tmax
                if not (b1 and b2):
                    continue
            X.append(x/2**32*360)
            Y.append(y/2**32*360)#2^32/360
    X = np.asarray(X)
    Y = np.asarray(Y)
    return X,Y
    #plt.plot(X,Y,'b.')

def compare_times(t1,t2,comp = '>'):
    if comp == '>':
        t1.i
        
def get_traj(filename,tmin=None,tmax=None):
    d = decode(filename)
    X,Y = convert_traj(d,tmin=tmin,tmax=tmax)
    return X,Y
