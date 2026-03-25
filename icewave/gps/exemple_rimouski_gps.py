

#non standard modules :
import gpxpy
try:
    import tilemapbase #unstable module
except:
    print("Cannot import module tilemapbase")
    
import glob
import os
import numpy as np
import pylab as plt
import socket
import platform
import pickle

#personal modules :
import stephane.rimouski.gps as gps
import stephane.rimouski.garmin as garmin
import stephane.display.graphes as graphes

#Global variables
global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()

def exemples():
    base = 'Banquise/Rimouski_2023/Data/GPS/'
    date = '20230313'

    # get a trajectory, and plot it
    filename = get_files(base,date,ext='fit')[0]
    Long,Lat = garmin.get_traj(filename)
    ax,figs = gps.display_traj(filename,Long,Lat,save=True)
    if 'thiou' not in osname:
        plt.show()
        
    # get the waypoints, and display them on a local map
    filename = get_files(base,date,ext='gpx')[0]
    gpx = gps.get_wpts(filename)
    ax,figs = gps.display_gpx(filename,date,gpx,save=True)
    if 'thiou' not in osname:
        plt.show()
#    savefolder = os.path.dirname(filename)+'/'
#    graphes.save_figs(figs,savedir=savefolder,prefix='stereo',suffix='labeled',frmt='pdf')

    #load a dictionnary of waypoints, plot it on top of a situation map.
    fileGPS = get_files(base,date,ext='pickle',prefix='Data_GPS')[0]
    f = open(fileGPS, 'rb')
    wpts = pickle.load(f)
    ax,figs = gps.display_dictwpts(fileGPS,date,wpts,save=True)
    if 'thiou' not in osname:
        plt.show()

def get_files(base,date,ext=None,prefix=''):
    serveurfolder = find_path(base)
    print(serveurfolder+date+'/'+prefix+'*.'+ext)
    filelist = glob.glob(serveurfolder+date+'/'+prefix+'*.'+ext)
    print(filelist)
    return filelist    

def find_path(base):
	print('OS type : '+str(osname))

	if 'Windows' in ostype:    
		serveurfolder = 'W:/'+base #fucking windows OS : beware of the mounting disk you used
    
	if 'Linux' in ostype:
		if 'thiou' in osname:
			serveurfolder = '/volume3/labshared2/'+base #praise UNIX system    		
		else:
			serveurfolder = '/media/turbots/DATA/thiou/labshared2/'+base #praise UNIX system
	if 'Darwin' in ostype:# or 'laita' in osname:    
		serveurfolder = '/Volumes/labshared2/'+base #praise UNIX system        

	return serveurfolder

def get_dict():
    pass

if __name__== '__main__':
    exemples()
