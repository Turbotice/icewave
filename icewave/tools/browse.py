
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
		serveurfolder = 'W:/'+base #fucking windows OS : beware of the mounting disk you used
    
	if 'Linux' in ostype:
		if 'thiou' in osname:
			serveurfolder = '/volume3/'+disk+'/'+base #praise UNIX system    		
		else:
			serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base #praise UNIX system
	if 'Darwin' in ostype:# or 'laita' in osname:    
		serveurfolder = '/Volumes/'+disk+'/'+base #praise UNIX system        

	return serveurfolder
	
