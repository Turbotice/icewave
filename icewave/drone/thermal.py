#from PIL import Image
#from PIL.ExifTags import TAGS

from PIL import Image, ExifTags
from exiftool import ExifToolHelper
import base64
import struct

import matplotlib.pyplot as plt
import numpy as np
import glob

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw

from skimage import measure
import scipy.ndimage as nd


import os

def position_exif(exif):
    lat = exif['GPSInfo'][2]#old version of the exif
    latitude = float(lat[0]+lat[1]/60+lat[2]/3600)

    lon = exif['GPSInfo'][4]
    longitude = -float(lon[0]+lon[1]/60+lon[2]/3600)
    return latitude,longitude

def analyse_IR(filename):
    # need to add metadata to retrieve temperature
    im = plt.imread(imagename)

    #compute the histogram
    n,xc = histogram(im,N=100,Max=400)


def regionprops(bw,Thres,minsize=50):
    binary = np.asarray(bw)>Thres
    
    struct = [[1,1,1],[1,1,1],[1,1,1]]
    labels = nd.label(binary, structure=struct)[0]
    #colormap of the fracture, color codes the time of opening (in rot max ?)
    props = measure.regionprops_table(labels,properties=['area'])#,'perimeter','axis_major_length','axis_minor_length','centroid','eccentricity','intensity_mean','intensity_max','intensity_median','intensity_std'])

    minsize = 50 #remove all areas smaller than 50 pixels
    labnums = np.where(props['area']>=minsize)[0]+1
    #print(labnums)
    [ny,nx] = binary.shape

    bigs_mask = np.asarray([[(labels[i,j] in labnums) for j in range(nx)] for i in range(ny)])

    bigs = binary*bigs_mask
    labels = nd.label(np.asarray(bigs), structure=struct)[0]
    props = measure.regionprops_table(labels, intensity_image=bw,properties=['area','perimeter','axis_major_length','axis_minor_length','centroid','eccentricity','intensity_mean','intensity_max','intensity_median','intensity_std'])

    return props,bigs


def histogram(im,N=100,Max=400):
    norm = np.linalg.norm(im,axis=2)
    [n,x] = np.histogram(np.ndarray.flatten(norm),bins=N,range=[0,Max])
    xc = (x[1:]+x[:-1])/2
    return n,xc


def get_exif(filename,display=False,movie=False):
    exif = {}
    with ExifToolHelper() as et:
        if movie:
            metadata = et.get_metadata(filename,params = '-ExtractEmbedded')
        else:
            metadata = et.get_metadata(filename,params = '-b')

        for d in metadata:
            for k, v in d.items():
                if k=='QuickTime:Text':
                    print(v)
                if display:
                    print(f"Dict: {k} = {v}")
                exif[k] = v
    return exif

def get_angles_distances(data):
    img = data['IR']
    H = data['H']
    alpha_0 = data['pitch']
    focale = data['focale']

    Xreal,Yreal = georectify_image(img,H,alpha_0,focale)

def compute_hist_raw(raw):
    cmin = np.min(np.ndarray.flatten(raw))
    cmax = np.max(np.ndarray.flatten(raw))
    [n,x] = np.histogram(np.ndarray.flatten(raw),range=[cmin-1/2,cmax+1/2],bins=int((cmax-cmin)+1))
    xc = (x[1:]+x[:-1])/2
    return n,xc

def compute_IR_avg(data):
    #compute some quantities from raw data
    if data is None or not('APP3:ThermalData' in data.keys()):
        print('Cannot load Raw thermal data, skip')
        return None
    raw = data['APP3:ThermalData']
    data['IR_moy'] = np.mean(raw)
    data['IR_std'] = np.std(raw)

    n,xc = compute_hist_raw(raw)
    data['IRhist_n']=n
    data['IRhist_xc']=xc

    return data



def load_image(filename):
    exif = get_exif(filename)
    try:
        data_IR = read_rawdata(exif)
    except:
        print('Cannot load Raw data, skip')
        return None
    #list the fields to retrieve
    keys = {'File:FileName',
            'EXIF:DateTimeOriginal',
            'EXIF:ApertureValue',
            'EXIF:DigitalZoomRatio',
            'EXIF:XPComment',
            'EXIF:GPSLatitudeRef',
            'EXIF:GPSLatitude',
            'EXIF:GPSLongitudeRef',
            'EXIF:GPSLongitude',
            'EXIF:GPSAltitudeRef',
            'EXIF:DateTimeOriginal',
            'EXIF:FocalLength',
            'EXIF:GPSAltitude',
            'EXIF:ThumbnailOffset',
            'EXIF:ThumbnailLength',
            'EXIF:XPKeywords',
            'XMP:AbsoluteAltitude',
            'XMP:RelativeAltitude',
            'XMP:GPSLatitude',
            'XMP:GPSLongitude',
            'XMP:GimbalRollDegree',
            'XMP:GimbalYawDegree',
            'XMP:GimbalPitchDegree',
            'XMP:FlightRollDegree',
            'XMP:FlightYawDegree',
            'XMP:FlightPitchDegree',
            'XMP:FlightXSpeed',
            'XMP:FlightYSpeed',
            'XMP:FlightZSpeed',
            'XMP:SensorTemperature',
            'Composite:GPSAltitude',
            'Composite:GPSLatitude',
            'Composite:GPSLongitude',
            'Composite:CircleOfConfusion',
            'Composite:FocalLength35efl',
            'Composite:GPSPosition',
            'Composite:HyperfocalDistance'}
    data = {}
    data['APP3:ThermalData']=data_IR
    for key in keys:
        data[key]=exif[key]
    return data

def extract_image(filename):
    im =  Image.open(filename)
    nx = im.height
    ny = im.width
    im = np.reshape(np.asarray(im.get_flattened_data()),(nx,ny,3))#.shape
    return im

def extract_exif_mp4(filename):
    import base64
    exif = th.get_exif(filename,display=False,movie=False)
    base64_str = exif['QuickTime:PreviewImage']
    if base64_str.startswith("base64:"):
        base64_str = base64_str[7:]  # Enleve "base64:"
    # Décoder en bytes
    binary_data = base64.b64decode(base64_str)
    uint16_array_np = np.frombuffer(binary_data, dtype=np.uint16)
    im = extract_image(BytesIO(uint16_array_np))
    return im

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

def compare_RGB_images(im1,im2,label1='',label2='',vmin=100,vmax=150):
    fig,axs = plt.subplots(figsize=(12,10),ncols=2,nrows=2)
    
    sc = axs[0,0].imshow(im1,vmin=vmin,vmax=vmax)
    plt.colorbar(sc,location='bottom',ax=axs[1,0])
    sc = axs[0,1].imshow(im2,vmin=vmin,vmax=vmax)
    plt.colorbar(sc,location='bottom',ax=axs[1,1])
    
    axs[0,0].set_title(label1, fontsize=20)
    axs[0,1].set_title(label2, fontsize=20)
    
    
    #img.get_flattened_data()
    colors = ['r','g','b']
    for i,color in enumerate(colors):
        out1 = axs[1,0].hist(np.ndarray.flatten(im1[...,i]),range=[0,2**8],bins=2**8,facecolor=color)
        out2 = axs[1,1].hist(np.ndarray.flatten(im2[...,i]),range=[0,2**8],bins=2**8,facecolor=color)
    
    graphes.legende('R,G,B','$n$','',ax=axs[1,0])
    figs = graphes.legende('R,G,B','$n$','',ax=axs[1,1])
    return figs


def compare_RGB_raw(im,imraw,label1='',label2='',vmin=100,vmax=150):
    fig,axs = plt.subplots(figsize=(12,10),ncols=2,nrows=2)
    
    sc = axs[0,0].imshow(im,vmin=vmin,vmax=vmax)
    plt.colorbar(sc,location='bottom',ax=axs[1,0])
    sc = axs[0,1].imshow(imraw)
    plt.colorbar(sc,location='bottom',ax=axs[1,1])
    
    axs[0,0].set_title(label1, fontsize=20)
    axs[0,1].set_title(label2, fontsize=20)
    
    
    cmin = np.min(np.ndarray.flatten(imraw))
    cmax = np.max(np.ndarray.flatten(imraw))
    #img.get_flattened_data()
    colors = ['r','g','b']
    for i,color in enumerate(colors):
        out1 = axs[1,0].hist(np.ndarray.flatten(im[...,i]),range=[0,2**8],bins=2**8,facecolor=color)
    out2 = axs[1,1].hist(np.ndarray.flatten(imraw),range=[cmin,cmax],bins=int((cmax-cmin)),facecolor=color)

    graphes.legende('R,G,B','$n$','',ax=axs[1,0])
    figs = graphes.legende('Uint16','$n$','',ax=axs[1,1])
    return figs

def exemple_2():
    folder = '/Users/stephane/Documents/BicWin2026/IR_Images_Buffer/'
#    folder = '/Users/stephane/Documents/BicWin2026/IR_data/IR_Images_Buffer/'
    filelist = glob.glob(folder+'*_T.JPG')
    print(filelist)

    datas = {}
    for filename in filelist:
        print(filename)
        name = os.path.basename(filename).split('.')[0]
        data = load_image(filename)
        if data is not None:
            datas[name] = data
    return datas

def gen_all():
    #print(read_rawdata())
    base = df.find_path(disk='Shack26')
    filelist = glob.glob(base+'*/Drones/Tourterelle/*/*_T.JPG')
    folders = set([filename.split('/')[-2] for filename in filelist])

    savefolder = base + 'IR_Data/extract_RJPG/'
    for folder in folders:
        datas = {}
        filelist = glob.glob(base+f'*/Drones/Tourterelle/{folder}/*_T.JPG')
        for filename in filelist:
            print(filename)
            name = os.path.basename(filename).split('.')[0]
            data = load_image(filename)
            data = compute_IR_avg(data)
            if data is None:
                continue
            data.pop('APP3:ThermalData', None)        
            datas[name] = data
    #rw.
        namef = os.path.basename(folder)
        filename = savefolder +f'{namef}_IR_datas_histograms.pkl'
        rw.write_pkl(filename,datas)
    return datas


def exemple_3():
    #adapt PATH to your local machine
    folder = '/Users/stephane/Documents/BicWin2026/Data_Exemples/IR_0211/'
    savefolder =  '/Users/stephane/Documents/BicWin2026/Data_Exemples/Results/'

    if not os.path.exists(savefolder):
        os.makedirs(savefolder)

        
    videoname = folder+'DJI_20260212054406_0055_T.MP4'
    rjpgname1 =  folder+'DJI_20260212054356_0054_T.JPG'
    rjpgname2 =  folder+'DJI_20260212054715_0056_T.JPG'

    extractimage1 = folder+'extract/frame_0001.PNG'
    extractimage2 = folder+'extract/frame_5419.PNG'


    imvid1 = extract_image(extractimage1)
    imvid2 = extract_image(extractimage2)

    imrjpg1 = extract_image(rjpgname1)
    imrjpg2 = extract_image(rjpgname2)

    #extract first image of MP4 from exif
    imvid1_exif = extract_exif_mp4(videoname)

    #extract raw data from exif of RJPG
    exif = th.get_exif(rjpgname1,display=False,movie=False)
    rawdata1 = th.read_rawdata(exif)

    exif = th.get_exif(rjpgname2,display=False,movie=False)
    rawdata2 = th.read_rawdata(exif)


    figs = compare_RGB_images(imvid1,imrjpg1,label1='From ffmpeg',label2='From _T.JPG')
    graphes.save_figs(figs,savedir=savefolder,prefix='MP4_RJPG_ffmpeg_RJPG',overwrite=True)
    

if __name__=='__main__':
    gen_all() 
