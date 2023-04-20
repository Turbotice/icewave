import glob
import pylab
import os

import imageio
import cv2

import numpy as np

import icewave.tools.browse as browse

def main():
    date = '20230310'
    datafolder = 'Banquise/Rimouski_2023/Data/drone/'+date+'/'
    datafolder = browse.find_path(datafolder,disk='labshared2')

    print(datafolder)
    file_mov = glob.glob(datafolder+'contexte/video/*.MOV')
    file_mp4 = glob.glob(datafolder+'contexte/video/*.MP4')
    	
    filelist = file_mov+file_mp4
    print(len(filelist))

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
