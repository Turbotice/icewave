#%%
import numpy as np
import pickle

import csv
from datetime import datetime
import pytz
import glob 
import scipy
import sys

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/')

from icewave.tools import weather

#%%

interp_H = weather.get_bathy_interpolator(disk='D:/copie_BicWin25_geophones',year='2025')


# %% importer les données gps des acquisitions 
# de géophones et en déduire la bathymétrie à ces endroits

#csv_file_coord = 'B:/General/Summary_geophone_lines/' + 'all_geophone_lines_coordinates.csv'

csv_file_coord = 'D:/copie_BicWin25_geophones/General/Summary_geophone_lines/' + 'all_geophone_lines_coordinates.csv'

with open(csv_file_coord, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data = []
    for lines in csvFile:
        if count==0:
            header = lines
        else:
            data.append(lines)
        count+=1

print(data)    

list_H = []
for i in range(len(data)):
    H = float(weather.get_bathymetry_GPS((float(data[i][3]), float(data[i][2])), interp_H))
    list_H.append(H)
    # format : latitude, longitude
print(list_H)
# %%
# à continuer... (avoir aussi les heures des acquisitions pour trouver les marées...)
#csv_file_times_UTC = 'B:/General/Summary_geophone_lines/' + 'all_geophone_lines_times_UTC.csv'
csv_file_times_UTC = 'D:/copie_BicWin25_geophones/General/Summary_geophone_lines/' + 'all_geophone_lines_times_UTC.csv'

with open(csv_file_times_UTC, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data1 = []
    for lines in csvFile:
        if count==0:
            header1 = lines
        else:
            data1.append(lines)
        count+=1
print(header1)

print(data1) 

acq_numbers = []
tide_height_array = []
yyyymmdd_array = []
UTC_datetime_array = []
for i in range(len(data1)):
    acq_numbers.append(int(data1[i][1]))
    UTC_datetime_str = data1[i][2] # on prend l'instant où les géophones ont été allumés
    print(UTC_datetime_str)
    UTC_datetime = datetime.fromisoformat(UTC_datetime_str.replace("Z", "+00:00"))
    #tide_height = weather.tide_from_datetime(UTC_datetime,disk = 'B:',year = '2025')
    tide_height = weather.tide_from_datetime(UTC_datetime,disk = 'D:/copie_BicWin25_geophones',year = '2025')
    print(tide_height)
    tide_height_array.append(tide_height)
    yyyymmdd_array.append(str(UTC_datetime.year).zfill(2)+str(UTC_datetime.month).zfill(2)+str(UTC_datetime.day).zfill(2))
    UTC_datetime_array.append(UTC_datetime)
#%% additionner pour avoir les profondeurs de chaque acquisition

water_depth_array = np.array(list_H) + np.array(tide_height_array)


# %% aller chercher dans data les coordonnées des géophones pour 
# chaque acquisition référencée dans data1

data2save = []

for i in range(len(data1)):
    f=0
    if (data1[i][0]==data[i][0])&(data1[i][1]==data[i][1]):
        data2save.append([yyyymmdd_array[i], acq_numbers[i], data[i][2], data[i][3], water_depth_array[i], data1[i][2]])
# %% save the file save2file

savefolder = 'D:/copie_BicWin25_geophones/a_copier_sur_Backup25/General/Summary_geophone_lines/water_level_acquisitions/'
file2save_path = savefolder + 'coord_waterlevel_allacq.csv'

header2save = ['date', 'acq', 'longitude', 'latitude', 'water level (m)', 'UTCtime_acqstart']

with open(file2save_path, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f,delimiter=',')

    # write the header
    writer.writerow(header2save)

    # write multiple rows
    writer.writerows(data2save)

# %%
