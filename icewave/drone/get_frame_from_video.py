# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 11:10:55 2025

@author: sebas

This script extract a frame of a video and save it in the same folder 

"""

import glob
# import pylab
import imageio
import cv2
import numpy as np
import os
import icewave.tools.datafolders as df


def main():
    # put your correct main_path
    main_path = 'K:/Share_hublot/Repository/'
    
    # replace these parameters depending on the video you want to extract a frame 
    date = '0306'
    drone_ID = 'bernache'
    exp_ID = '12-doc_vuedensemble'
    video_suffixe = '*0198_D.MP4'
    idx_frame = 60
    
    path = f'{main_path}{date}/Drones/{drone_ID}/{exp_ID}/{video_suffixe}'
    save_frame_from_video(path, idx_frame)
    print('DONE.')


def save_frame_from_video(path,idx_frame):
    """ Save frame number idx_frame from a video, saved under path """
    file2load = glob.glob(path)[0]
    print(file2load)
    
    vid = imageio.get_reader(file2load,  'ffmpeg')
    image = vid.get_data(idx_frame)
    
    frame_name = file2load.replace('.MP4',f'_frame_{idx_frame}.tiff')
    cv2.imwrite(frame_name, cv2.cvtColor(image, cv2.COLOR_RGB2BGR))


if __name__=='__main__':
    
    main()
