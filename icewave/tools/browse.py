
import glob
import os
#Global variables0

import platform
import socket

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
    
    
