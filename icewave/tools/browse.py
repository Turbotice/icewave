
import glob
import os
#Global variables0

import platform
import socket

global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()


def create_folder(folder):
    #also adjust automatically the rights on the created folder ?
    if not os.path.exists(folder):
        print('Warning, '+folder+' does not exist. Do you want to create it ?')
        os.makedirs(folder)

def find_path(base,disk='labshared2'):
    print('OS type : '+str(osname))

    if 'Windows' in ostype:    
        if disk=='labshared2':
            serveurfolder = 'W:/'+base #fucking windows OS : beware of the mounting disk you used
        if disk=='labshared1':
            serveurfolder = 'Y:/'+base #fucking windows OS : beware of the mounting disk you used
        if disk=='storageshared':
            serveurfolder = 'X:/'+base #fucking windows OS : beware of the mounting disk you used
        if disk=='homes':
            serveurfolder = 'Z:/'+base #fucking windows OS : beware of the mounting disk you used

    if 'Linux' in ostype:
        if 'adour' in osname:
            serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base
        elif 'spi201711-Latitude-5480' in osname:
            serveurfolder = '/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/'+disk+'/'+base
        elif 'thiou' in osname:
            serveurfolder = '/volume3/'+disk+'/'+base #praise UNIX system    		
        elif 'oural' in osname:
            serveurfolder = '/media/turbots/'+disk+'/'+base

        else:
            serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base #praise UNIX system
			
    if 'Darwin' in ostype:# or 'laita' in osname:    
        serveurfolder = '/Volumes/'+disk+'/'+base #praise UNIX system        

    return serveurfolder	

def folders(date,folder='Banquise/Rimouski_2023/Data/Geophones/'):
    savefolder = '/Users/stephane/Documents/Programming/Python/Notebooks/Geophone_Rimouski/'
    #    folder = '/Users/stephane/Documents/Rimouski 2023/Data/Geophones/'+date+'//'
    #    savefolder = '/Users/stephane/Documents/Rimouski 2023/Data/Geophones/20230313/Resultats/'
    folder = folder+date+'/'
    
    if 'Linux' in ostype:
        if 'adour' in osname:
            base = '/media/turbots/DATA/thiou/labshared2/'
        if 'spi201711-Latitude-5480' in osname:
            base = '/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/labshared2/'
            
    if 'Darwin' in ostype:# or 'laita' in osname:    
        base = '/Volumes/labshared2/' 

    savefolder = base + folder + 'Results_Sarah/'

    if not os.path.exists(savefolder):
        os.makedirs(savefolder)

    return base,folder,savefolder
    
    

