# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 16:48:05 2024

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
    datafolder = '/media/turbots/DATA/thiou/labshared2/SagWin2024/Data/Drone/2023/0927/0001/'
    
    # datafolder = browse.find_path(datafolder)
    print(datafolder)
    file_mp4 = glob.glob(datafolder+'*.MP4')
    
    filelist = file_mp4
    print(len(filelist))
    print(filelist)
    return filelist


# for filename in filelist[2]:
#     print(os.path.basename(filename))
def later(filelist,datafolder):
    filename = filelist[0]
    print(filename)
    vid = imageio.get_reader(filename,  'ffmpeg')
    nb_frames = vid.count_frames() # number of frames
    nums = [p for p in range(0,nb_frames)]
    #directory = os.path.dirname(filename)
    directory = datafolder
    print(directory)
    
    savefolder = directory +'/'+os.path.basename(filename).split('.')[0] +'/'
    print(savefolder)
    if not os.path.isdir(savefolder):
        os.makedirs(savefolder)
            
    for num in nums:
        print(num)
        image = vid.get_data(num)
        # print(image.shape)
        
        nb_max_digit = len(str(nums[-1]))
        diff = nb_max_digit - len(str(num))
        string_index = ''
        for i in range(diff):
            string_index = '0' + string_index
        string_index = string_index + str(num)
        fname = savefolder + '/im_'+ string_index +'.tiff'
        cv2.imwrite(fname,image)
        
        #     fig = pylab.figure()
        #     fig.suptitle('image #{}'.format(num), fontsize=20)
        #     pylab.imshow(image)
        # pylab.show()


if __name__=='__main__':
    filelist = main()
    datafolder = '/media/turbots/DATA/thiou/labshared2/SagWin2024/Data/Drone/2023/0927/0001/'
    vid = later(filelist,datafolder)
