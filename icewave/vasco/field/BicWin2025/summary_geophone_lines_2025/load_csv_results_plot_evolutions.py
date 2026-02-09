#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
from datetime import datetime
import matplotlib.dates as mdates

#%%
disk = 'B:'
csvfilepath = f"{disk}/General/Summary_geophone_lines/results_acquisitions.csv"


with open(csvfilepath, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data_csv = []
    for lines in csvFile:
        if count==0:
            header = lines
        elif count>=1:
            data_csv.append(lines)
        count+=1

data_csv = np.array(data_csv)

print(header)

#%%
def convert_array1d_to_floatarray(arr):
    arr_new = []
    for i in range(len(arr)):
        if arr[i]=='':
            arr_new.append(np.nan)
        else:
            arr_new.append(float(arr[i]))
    return np.array(arr_new)


def convert_array1d_to_datetimearray(arr):
    arr_new = []
    for i in range(len(arr)):
        dt = datetime.strptime(arr[i], "%Y-%m-%dT%H:%M:%S.%fZ")
        arr_new.append(dt)
    return np.array(arr_new)   

# %%
%matplotlib qt

plt.figure()
plt.plot(convert_array1d_to_datetimearray(data_csv[:,5]), convert_array1d_to_floatarray(data_csv[:,7]), '.b')
plt.ylim(0, 7e9)
plt.show()
# %% fonction permettant de specifier des manips de geophones dans une certaine zone : exemple baie du ha!ha!


def points_in_bbox(lon, lat, west, east, north, south):
    """
    Vérifie quels points sont à l'intérieur d'une boîte géographique.

    Paramètres
    ----------
    lon : array-like
        Longitudes des points (1D)
    lat : array-like
        Latitudes des points (1D)
    west : float
        Longitude minimale (limite ouest)
    east : float
        Longitude maximale (limite est)
    north : float
        Latitude maximale (limite nord)
    south : float
        Latitude minimale (limite sud)

    Retour
    ------
    inside : np.ndarray (bool)
        Tableau booléen True/False indiquant si chaque point est dans la boîte.
    """

    lon = np.asarray(lon)
    lat = np.asarray(lat)

    inside = (lon >= west) & (lon <= east) & (lat >= south) & (lat <= north)
    return inside

west_haha = -68.85366812475961
east_haha = -68.80077777777778
north_haha = 48.3564444
south_haha = 48.33393038416758


#%%
longitude_arr = convert_array1d_to_floatarray(data_csv[:,2])
latitude_arr = convert_array1d_to_floatarray(data_csv[:,3])
datetime_arr = convert_array1d_to_datetimearray(data_csv[:,5])
h_arr = convert_array1d_to_floatarray(data_csv[:,6])
E_arr = convert_array1d_to_floatarray(data_csv[:,7])
nu_arr = convert_array1d_to_floatarray(data_csv[:,8])
rho_arr = convert_array1d_to_floatarray(data_csv[:,9])


#indices_located = points_in_bbox(longitude_arr, latitude_arr, west_haha, east_haha, north_haha, south_haha)

west = -68.818
east = -68.81326448878322
north = 48.34954246746782
south = 48.34687414439287

indices_located = points_in_bbox(longitude_arr, latitude_arr, west, east, north, south)

fig, axs = plt.subplots(2, 2, figsize=(10, 6))

axs[0,0].plot(datetime_arr[indices_located], h_arr[indices_located], 'o')
axs[0,1].plot(datetime_arr[indices_located], E_arr[indices_located], 'o')
axs[1,0].plot(datetime_arr[indices_located], nu_arr[indices_located], 'o')
axs[1,1].plot(datetime_arr[indices_located], rho_arr[indices_located], 'o')


axs[0,0].set_title('thickness (m)')
axs[0,1].set_title('Young modulus (Pa)')
axs[1,0].set_title('Poisson ratio')
axs[1,1].set_title('Density (kg/m$^3$)')


# Boucle de formatage des 4 axes
for i in range(2):
    for j in range(2):
        ax = axs[i, j]
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m %H:%M'))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        # rotation et alignement pour chaque subplot
        for label in ax.get_xticklabels():
            label.set_rotation(30)
            label.set_ha('right')

plt.tight_layout()
plt.show()


# %%
