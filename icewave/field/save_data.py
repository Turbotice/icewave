
import icewave.gps.garmin as gar
import icewave.pyphone_v2.ls_phone as phone


from pprint import pprint
#save locally all data of the day, store it directly in the tree of folders
#

def download():
    gar.download()
    phone.download()

if __name__=='__main__':
    download()
