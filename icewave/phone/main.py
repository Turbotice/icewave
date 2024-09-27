
import glob
import os
#Global variables0

import platform
import socket
import zipfile as zip
import rw_pyphone as rw

global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()


def main():
    testfolder = 'Telephones/Soufflerie_dec23/131223/Telephones/121223_4_U400cms/'

    basefolder = 'Telephones/Soufflerie_dec23/131223/Telephones/
#    rw.extract_all(testfolder)
    folder = rw.find_path(testfolder)
    csvlist = glob.glob(folder+'*/*.csv')

#    filelist = rw.get_ziplist(testfolder)
    #print(filelist)
    for filename in csvlist:
        print(rw.get_phone_num(filename), rw.get_filename_key(filename))

    return csvlist

if __name__ =='__main__':
    main()
