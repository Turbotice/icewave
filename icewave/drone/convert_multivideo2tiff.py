# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:32:58 2024

@author: sebas
"""
import glob
# import pylab
import imageio
import cv2
import numpy as np
import os
import icewave.tools.datafolders as df

#%%

def main(base,date,drone_ID,exp_ID):
    # base = df.find_path('Hublot24')
    datafolder = base + f'{date}/Drones/{drone_ID}/{exp_ID}/'

    # datafolder = browse.find_path(datafolder)
    print(datafolder)
    file_mp4 = glob.glob(datafolder+'*.MP4')
    
    filelist = file_mp4
    filelist = sorted(filelist)
    print(len(filelist))
    print(filelist)
    return filelist

def find_number_frames(filelist):
    
    N = 0
    frames_array = np.zeros(len(filelist))
    
    for i in range(len(filelist)) :
        filename = filelist[i]
        print(filename)
        vid = imageio.get_reader(filename,  'ffmpeg')
        nb_frames = vid.count_frames() # number of frames in the video
        if i == 0 :
            frames_array[i] = nb_frames
        else:
            frames_array[i] = int(frames_array[i-1] + nb_frames)              
    N = int(frames_array[-1]) 
    print(N)
    frames_array = frames_array.astype(int)
    print(frames_array)
    return N, frames_array
    

# for filename in filelist[2]:
#     print(os.path.basename(filename))
def later(filelist,directory_save,savename):
    
    [N,frames_array] = find_number_frames(filelist) # total number of frames
    nb_max_digit = len(str(N)) 
    
    savefolder = directory_save +'/'+os.path.basename(savename).split('.')[0] +'/'
    print(savefolder)
        
    if not os.path.isdir(savefolder):
        os.makedirs(savefolder)

    for i in range(len(filelist)):

        filename = filelist[i]
        print(filename)
        vid = imageio.get_reader(filename,  'ffmpeg')
         
        if i == 0 :
            nums = [p for p in range(0, frames_array[0])]
        else :
            nums = [p for p in range(frames_array[i-1],frames_array[i])]
            
        print(nums[-1])
    #directory = os.path.dirname(filename)
       
        for num in nums:
            print(num)
            idx = num - nums[0]
            image = vid.get_data(idx)
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

    print('DONE.')    

if __name__=='__main__':

    # base = df.find_path('Hublot24')

    base = 'F:/Rimouski_2025/Data/'
    date = '0210'
    drone_ID = 'bernache'
    exp_ID = '11-situation_map_geophone_003/movie_1'
    filelist = main(base,date,drone_ID,exp_ID)
    directory_save = base[:-5] + f'/PIV_images/{date}/Drones/{drone_ID}/'
    vid = later(filelist,directory_save,exp_ID)
