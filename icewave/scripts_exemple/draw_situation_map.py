
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

import icewave.drone.drone_projection as proj
import icewave.display.maps as maps
import icewave.display.graphes as graphes

global date
date = '0223'
ti = '17:50:00' #initial time in UTC time
tf = '18:40:00' #final time in UTC time


base = df.find_path()
print(base+date)
savefolder = base+'Summary/Results/'

filename =  base+f'{date}/Summary/records_{date}.pkl'
records = rw_data.load_pkl(filename)
print(records.keys())

#tifflist = glob.glob(base+date+'/Drones/Bernache/*/*.tiff')
#print(tifflist)
drone = 'Bernache'
filename =  base+f'{date}/drones/{drone}/flightrecords/Flightrecord_dict.pkl'
flight = rw_data.load_pkl(filename)
print(flight.keys())

keydrone = '20-waves_013'
record = records['drones'][drone][keydrone][0]
flight_p = drone.cut_flightrecord(record,flight)

im = proj.get_exemple_image(record)

Lats,Lons = proj.project_image(record,flight_p,im,0)

# draw map
fig,ax = plt.subplots(figsize=(15,4))
ax.pcolormesh(Lons, Lats, im)

#ax.plot(vertices_real[:,0,:],vertices_real[:,0,:],'go')       
measures = maps.get_measures(records,ti,tf)
figs = maps.display_map(measures,remote=True,ax=ax)#w=10,h=6,
figs = graphes.legende('Longitude','Latitude',f'{date}, {ti} to {tf} UTC')

ax.set_ylim([maps.coord2angle(48,19,45),maps.coord2angle(48,19,48)])
#graphes.save_figs(figs,savedir=savefolder,prefix=f'map_{date}_',overwrite=True,frmt='png')
#graphes.save_figs(figs,savedir=savefolder_local,prefix=f'map_{date}_',overwrite=True,frmt='png')
