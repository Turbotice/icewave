import glob
import pylab
import os

import imageio
import cv2

import numpy as np

import icewave.tools.browse as browse
import icewave.tools.datafolders as df

import subprocess

def main():
    base = df.find_path(year='2026')
    savefolder = base+'Buffer/IR_images/'
   
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    print(base)
    filelist = glob.glob(base+'*/Drones/Tourterelle/*/*_T.JPG')
    print(len(filelist))
    folder = []
    for filename in filelist:
        f = os.path.dirname(filename)
        if not f in folder:
            folder.append(f)
            print(filename)
            dest = savefolder + os.path.basename(filename)
            subprocess.run(['cp',filename,dest])


def other():
    for filename in file_mov+file_mp4:
        print(os.path.basename(filename))

    for filename in filelist:
        vid = imageio.get_reader(filename,  'ffmpeg')
        
        nums = range(1000)#np.arange([2**p for p in range(8)]
        directory = os.path.dirname(filename)
        print(directory)
        
        savefolder = directory +'/'+os.path.basename(filename).split('.')[0] + '_images'
        if not os.path.isdir(savefolder):
            os.makedirs(savefolder)
            
        
            
        for num in nums:
            fname = savefolder + '/im_'+str(num)+'.tiff'
            if os.path.exists(fname):#skip if the image already exists
            	continue

            try:
            	image = vid.get_data(num)
            except:
            	print('image number out of range')
            	break
            if np.mod(num,10)==0:
                print(num)
            

            cv2.imwrite(fname,image)
            
#        fig = pylab.figure()
#        fig.suptitle('image #{}'.format(num), fontsize=20)
#        pylab.imshow(image)
#    pylab.show()


if __name__=='__main__':
    main()
