#from PIL import Image
#from PIL.ExifTags import TAGS

from PIL import Image, ExifTags
from exiftool import ExifToolHelper
import base64

import matplotlib.pyplot as plt
import numpy as np


def position_exif(exif):
    lat = exif['GPSInfo'][2]
    latitude = float(lat[0]+lat[1]/60+lat[2]/3600)

    lon = exif['GPSInfo'][4]
    longitude = -float(lon[0]+lon[1]/60+lon[2]/3600)
    return latitude,longitude

def analyse_IR(filename):
    # need to add metadata to retrieve temperature
    im = plt.imread(imagename)

    #compute the histogram
    n,xc = histogram(im,N=100,Max=400)    


def histogram(im,N=100,Max=400):    
    norm = np.linalg.norm(im,axis=2)
    [n,x] = np.histogram(np.ndarray.flatten(norm),bins=N,range=[0,Max])
    xc = (x[1:]+x[:-1])/2
    return n,xc


def get_exif(filename):
    exif = {}
    with ExifToolHelper() as et:
        for d in et.get_metadata(filename,params = '-b'):
            for k, v in d.items():
                #print(f"Dict: {k} = {v}")
                exif[k] = v            
    return exif

def read_rawdata(exif):
    # Ta chaîne Base64
    base64_str = exif['APP3:ThermalData']
# Extraire la partie après "base64:"
    if base64_str.startswith("base64:"):
        base64_str = base64_str[7:]  # Enleve "base64:"
    # Décoder en bytes
    binary_data = base64.b64decode(base64_str)
    #print(f"Taille des données binaires : {len(binary_data)} octets")

    # Convertir en numpy array
    uint16_array_np = np.frombuffer(binary_data, dtype=np.uint16)

    #print(len(uint16_array_np))

    # Récupère les dimensions réelles du capteur
    nx = int(exif['File:ImageHeight']/2)
    ny = int(exif['File:ImageWidth']/2)

    # Conversion en image
    data = np.reshape(uint16_array_np,(nx,ny))
    return data


def exemple():
    folder = '/Users/stephane/Documents/BicWin2026/Data_Exemples/IR_0211/'
    # path to the image or video
    imagename = folder + "DJI_20260212054323_0043_T.jpg"

    # read the image data using PIL
    im = plt.imread(imagename)
    image = Image.open(imagename)


    #plt.imshow(im)
    #plt.show()

    # extract EXIF data
    exifdata = image.getexif()

    # iterating over all EXIF data fields
    for tag_id in exifdata:
        # get the tag name, instead of human unreadable tag id
        tag = TAGS.get(tag_id, tag_id)
        data = exifdata.get(tag_id)#.decode("utf-16")
        print(f"{tag:25}: {data}")  

