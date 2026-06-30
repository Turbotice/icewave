#from PIL import Image
#from PIL.ExifTags import TAGS

from PIL import Image, ExifTags


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
    image = Image.open(filename)
    exif = { ExifTags.TAGS[k]: v for k, v in image._getexif().items() if k in ExifTags.TAGS }
    return exif


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

