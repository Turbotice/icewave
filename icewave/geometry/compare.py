import numpy as np
from pprint import pprint
import glob


import icewave.geometry.display as display
import icewave.geometry.tables as tables
import icewave.geometry.save as save

global norme
date = '0223'
import icewave.geometry.experiment_2024_0223 as experiment

norme = display.read_norme()
global path_GPS
path_GPS = f"/Users/stephane/Documents/git/Data_local/{date}/GPS/"

def load():
    #load GPX table
    import gpxpy
    filelist = glob.glob(path_GPS+'*.gpx')
    filegpx = filelist[0]
    gpx_file = open(filegpx, 'r')
    gpx = gpxpy.parse(gpx_file)

    #load GPS table
    table = tables.read_table(path_GPS)
    ftab = np.asarray(table)
    #print(table[:,0])
    
    #load geometry table
    figs,geom = experiment.Sag24(display=False)
    fgeom = np.asarray(geom)
    
    print(fgeom[:,0])

    geom[0]+=['GPS_tag_#','GPS_time']
    for wpt in gpx.waypoints:
        num = str(int(wpt.name.split('Sag24')[1]))
        indices = np.where(num==ftab[:,0])[0]
        if len(indices)>0:
            tags = [table[i][1] for i in indices]
            print(num,tags)
            for tag in tags:
                indg = np.where(tag==fgeom[:,0])[0]
                if len(indg)>0:
                    j = indg[0] #geometry should correspond to only GPS point
                    print(num,indg,tag,geom[j][1:])
                    geom[j]+=[wpt.name]
                    s = wpt.time.strftime("%Y%m%dT%H%M%S")
                    print(s)
                    geom[j]+=[s]
                    geom[j]
                else:
                    pass
                    #print("Waypoint not found")

    ncol = len(geom[0])
    for i,g in enumerate(geom):
        if len(g)<ncol:
            geom[i]+=[' ' for i in range(ncol-len(g))]
    pprint(geom)

    filename = path_GPS+f'Geom_table_{date}_sync.txt'
    save.save(geom,filename)
    
if __name__ == '__main__':
    pass
