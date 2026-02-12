import time

import subprocess
from multiprocessing import Process


from pprint import pprint

import icewave.pyphone_v2.s_acqui as s_acqui
import icewave.pyphone_v2.run_phyphox as run_phyphox
import icewave.pyphone_v2.connect as connect


def folders():
    folder = 'storage/self/primary/Download/'
    foldersave = '/home/turbots/Documents/Bicwin2024/git/icewave/icewave/pyphone_v2/Params/'
    return folder,foldersave

def test1(num):
    folder,foldersave = folders()
    ls = ['adb','-s','192.168.0.'+str(num+100),'shell','ls',folder]
    output = subprocess.run(ls,stdout=subprocess.PIPE)
    pprint(str(output.stdout).split('\\n'))

    name = 'accelero_gyro_magneto_gps.phyphox'
    filename = folder+name

    pull = ['adb','-s','192.168.0.'+str(num+100),'pull',filename,foldersave+name]
    print(pull)
    output = subprocess.run(pull,stdout=subprocess.PIPE)

def test2(num):
    folder,foldersave = folders()
    name = 'accelero_gyro_magneto_gps.phyphox'
    push = ['adb','-s','192.168.0.'+str(num+100),'push',foldersave+name,folder+name]
    print(push)
    output = subprocess.run(push,stdout=subprocess.PIPE)
    print(output.out)
def test3(num):
    folder,foldersave = folders()
    name = 'accelero_gyro_magneto_gps.phyphox'
    update = ['adb','-s','192.168.0.'+str(num+100),'shell','am','start','--user','0','-a','android.intent.action.VIEW','-d','file:///'+folder+name,'-t','image/jpg','-p','de.rwth_aachen.phyphox']
    print(update)
    output = subprocess.run(update,stdout=subprocess.PIPE)

def test_50Hz(num,Tacq):

#    ip = '192.168.0.'+str(num+100)
    Tacq = 3600*2

    adrlist,phonelist = connect.connect()
    iplist = connect.get_adresslist(phonelist)

    tstamp = int(time.time())
    keywords = {'T': Tacq,'folder':'PhyPhox_'+str(tstamp)}
    p1 = Process(target=run_phyphox.run_serie, args=(iplist,),kwargs=keywords)
    p1.start()

    t0 = time.time()
    t1 = time.time()

    c=0
    T=60
    while (t1-t0)<Tacq:
        t1 = time.time()

        if (t1-t0)>c*T:
            c+=1
            out = run_phyphox.run_get(iplist,'all')

            pprint()
            pprint(int(t1-t0))
            pprint(out)        #for r in rlist['buffers'][1:4]:
        #s=s+r['name']+"&&"
        #s=s[:-2]
        time.sleep(T)

if __name__=='__main__':

    #test1(4)
#    test2(4)
#    test3(4)
    test_50Hz(1)
