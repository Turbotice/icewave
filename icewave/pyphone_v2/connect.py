

import subprocess
import time


global path
path = '/home/turbots/Documents/Bicwin2024/git/icewave/icewave/pyphone_v2/adb_usb_liste.txt'
global local_folder
localfolder = '/Users/stephane/Documents/PMMH/Telephones/test/'


def ipbase(network=0):
    return f'192.168.{network}.'

def basefolder():
    return 'Data_test/'

def get_adress(phone):
    return ipbase()+str(100+phone)

def get_adresslist(phonelist):
    adresslist = []
    for phone in phonelist:
        adresslist.append(get_adress(phone))
    return adresslist

def get_connected():
    idlist = connected()
    d = {}
    for id in idlist:
        phone = identify(id)
        if phone>=0:
            d[phone]={}
            d[phone]['id']=id
            print(phone,id)
    return d

def connected():
    c = subprocess.run(['adb','devices'],text=True,capture_output = True,shell=False)
    s = c.stdout.split('\n')

    idlist = []
    for line in s[1:]:
        if len(line)>0:
            id = line.split('\t')[0]
            idlist.append(id)
    return idlist

def identify(id):
    with open(path) as f:
        while True:
            line = f.readline()
            if len(line)<2:
                break
            num,pid = line.split('\n')[0].split('\t')
            #print(num,id)
            if pid==id:
                return int(num)
    print('id unknown, add it manually to adb_usb_list')
    return -1
    
def connect():
    #first test
    #c = subprocess.run(['ping',ip],text=True)
    folder = '/storage/self/primary/'

    c = subprocess.run(['adb','devices'],text=True,capture_output = True,shell=False)
    s = c.stdout.split('\n')
    adrlist,iplist,phonelist=[],[],[]

    for line in s[1:]:
        if line[:3]=='192':
            adr = line.split('\t')[0]
            
            print(adr)
            adrlist.append(adr)
            iplist.append(adr.split(':')[0])
            phone = int(adr.split(':')[0][-2:])
            
            if phone>=0 and phone<70:
	            phonelist.append(phone)
           
    print(phonelist)
    print("Number of phone connected :"+str(len(phonelist)))

    return adrlist,phonelist
    
def unlock(phone):
    address = get_adress(phone)

    c=subprocess.run(['adb','-s',address,'shell','input','keyevent','KEYCODE_WAKEUP'],text=True,capture_output=True)
    c=subprocess.run(['adb','-s',address,'shell','input','swipe','200','900','200','300'],text=True,capture_output=True)
    c=subprocess.run(['adb','-s',address,'shell','input','text','01012000'],text=True,capture_output=True)
    c=subprocess.run(['adb','-s',address,'shell','input','keyevent','66'],text=True,capture_output=True)
    
#adb -s ${adress} shell input swipe 200 900 200 300
#adb -s ${adress} shell input text 01012000
#adb -s ${adress} shell input keyevent 66
#done

def launch_app(phone,name):
    unlock(phone)
    time.sleep(0.2)
    address = get_adress(phone)
    
    app_table = {'phyphox':'de.rwth_aachen.phyphox/.ExperimentList',\
    		 'camera':'com.android.camera/.Camera'}
    
    if name in app_table.keys():
    	app = app_table[name]
    	print("Start "+str(name))
    	c=subprocess.run(['adb','-s',address,'shell','am','start-activity','-n',app],\
                          text=True,capture_output=True)
    else:
    	print("Activity not in list ! Please add it to connect_phone.app_table")
                              
def launch_camera(phone):
    unlock(phone)
    time.sleep(0.2)
    address = get_adress(phone)
    
    print("Start Camera")
    c=subprocess.run(['adb','-s',address,'shell','am','start-activity','-n','com.android.camera/.Camera'],\
                          text=True,capture_output=True)           
         
def take_picture(phone):
    address = get_adress(phone)
    c=subprocess.run(['adb','-s',address,'shell','input','keyevent','KEYCODE_CAMERA'],\
                          text=True,capture_output=True)
    
def launch_phyphox(phone):
    unlock(phone)
    time.sleep(0.2)
    address = get_adress(phone)
    # to go back to initial screen : 
    #adb -s 192.168.0.100 shell input keyevent 3    
    c=subprocess.run(['adb','-s',address,'shell','am','force-stop','-n','de.rwth_aachen.phyphox/.ExperimentList'],\
                          text=True,capture_output=True)
    time.sleep(0.2)
    c=subprocess.run(['adb','-s',address,'shell','am','start','-n','de.rwth_aachen.phyphox/.ExperimentList'],\
                          text=True,capture_output=True)
    #adb -s ${adress} shell am force-stop -n de.rwth_aachen.phyphox/.ExperimentList

def test_connect():
    c=subprocess.run(['adb','devices'],text=True,capture_output=True)
    log = c.stdout#.split('\n')
        #self.pt('toto')
    print(log)


if __name__=='__main__':
    launch_phyphox(0)
#    test_connect()
