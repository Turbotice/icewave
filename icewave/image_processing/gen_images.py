import glob
import pylab
import imageio


import icewave.tools.browse as browse


def main():
    datafolder = 'labshared2/Banquise/Rimouski_2023/Data/drone/'
    datafolder = browse.find_path(datafolder)

    file_mov = glob.glob(datafolder+'*/contexte/video/*.MOV')
    file_mp4 = glob.glob(datafolder+'*/contexte/video/*.MP4')

    filelist = file_mov+file_mp4
    print(len(filelist))

def later():
    for filename in file_mov+file_mp4:
        print(os.path.basename(filename))


    for filename in filelist:
        vid = imageio.get_reader(filename,  'ffmpeg')
        nums = [10]
        directory = os.path.dirname(filename)
        print(directory)
    
        for num in nums[:1000]:
            image = vid.get_data(num)
            print(image.shape)
        
#        fig = pylab.figure()
#        fig.suptitle('image #{}'.format(num), fontsize=20)
#        pylab.imshow(image)
#    pylab.show()


if __name__=='__main__':
    main()
