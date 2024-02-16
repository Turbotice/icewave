# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:32:58 2024

@author: sebas
"""
import glob
import pylab
import imageio
import cv2

import os
#import icewave.tools.browse as browse

#%%

def main():
    datafolder = '//192.168.1.70/Share/Data/0211/Drones/bernache/stereo_001/relevant/'

    # datafolder = browse.find_path(datafolder)
    print(datafolder)
    file_mp4 = glob.glob(datafolder+'*.MP4')
    
    filelist = file_mp4
    print(len(filelist))
    print(filelist)
    return filelist

def find_number_frames(filelist):
    
    
    N = 0
    for filename in filelist :
        print(filename)
        vid = imageio.get_reader(filename,  'ffmpeg')
        nb_frames = vid.count_frames() # number of frames in the video
        N += nb_frames 
        print(N)
    return N
    
    
# for filename in filelist[2]:
#     print(os.path.basename(filename))
def later(filelist,directory_save,savename):
    
    N = find_number_frames(filelist) # total number of frames
    nb_max_digit = len(str(N)) 
    
    savefolder = directory_save +'/'+os.path.basename(savename).split('.')[0] +'/'
    print(savefolder)
        
    if not os.path.isdir(savefolder):
        os.makedirs(savefolder)
    
    for i in range(len(filelist)):

        filename = filelist[i]
        print(filename)
        vid = imageio.get_reader(filename,  'ffmpeg')
    
        nb_frames = vid.count_frames() # number of frames in the video

        nums = [p for p in range(i, i + nb_frames)]
    #directory = os.path.dirname(filename)
       
        for num in nums:
            print(num)
            image = vid.get_data(num)
            # print(image.shape)
            
            
            diff = nb_max_digit - len(str(num))
            string_index = ''
            for i in range(diff):
                string_index = '0' + string_index
            string_index = string_index + str(num)
            fname = savefolder + '/im_'+ string_index +'.tiff'
            cv2.imwrite(fname, cv2.cvtColor(image, cv2.COLOR_RGB2BGR))
            #     fig = pylab.figure()
            #     fig.suptitle('image #{}'.format(num), fontsize=20)
            #     pylab.imshow(image)
            # pylab.show()


if __name__=='__main__':
    filelist = main()
    directory_save = '//192.168.1.70/Share/Data/0211/Drones/bernache/stereo_001/relevant/'
    savename = 'full_movie'
    vid = later(filelist,directory_save,savename)
