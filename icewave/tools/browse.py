
import glob
import os
#Global variables0
global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()


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
	
