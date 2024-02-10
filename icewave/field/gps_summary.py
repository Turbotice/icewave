
import icewave.gps.garmin as gar
import icewave.gps.gps as gps
import icewave.pyphone_v2.ls_phone as phone
import icewave.tools.datafolders as df

import glob
import os
from pprint import pprint
#save locally all data of the day, store it directly in the tree of folders
#
import stephane.display.graphes as graphes

def display_trajectories(date=''):
    if date == '':
        date = df.get_current_date()
    year,day = date.split('_')
    filelist = glob.glob('Bicwin2024/Data/'+year+'/'+day+'/GPS/*.fit')
    pprint('number of tracks : '+str(len(filelist)))
    
    for i,filename in enumerate(filelist):
        Long,Lat = gar.get_traj(filename)
        title = gps.title_std(filename,Long,Lat)
        print(title)
        
        ax,t,figs = gps.map_traj(Long,Lat,scale=1,save=False,title=title)
    #    ax,t,figs = gps.map_traj(Long,Lat,scale=1,save=False,title=title)
        savefolder = os.path.dirname(filename)+'/'
        graphes.save_figs(figs,savedir=savefolder)


if __name__=='__main__':

    display_trajectories('2024_0210')
