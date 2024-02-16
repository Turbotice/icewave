

import time
from multiprocessing import Process

import icewave.pyphone_v2.connect as connect
import icewave.pyphone_v2.run_phyphox as run_phyphox


import argparse

def run(Tacq):
    adrlist,phonelist = connect.connect()  
    base = connect.ipbase()
    iplist = [base + str(100+phone) for phone in phonelist]

    print('Run Phyphox on '+str(phonelist))
    tstamp = int(time.time())
    keywords = {'T': Tacq,'folder':'PhyPhox_'+str(tstamp)}
    
    t1 = Process(target=run_phyphox.run_serie, args=(iplist,),kwargs=keywords)
    t1.start()
    #battery_T.acquire(phone_dict,Tmax=100)
        #checked = Qt.CheckState(self.gridlayout.itemAt(i)) == Qt.CheckState.Checked
        #print(phone,checked)

if __name__ =='__main__':
    parser = argparse.ArgumentParser(description="Start Acquisition")
    parser.add_argument('-T',dest='Tacq',default=10,type=int,help='Temps d''acquisition ')
    args = parser.parse_args()

    if args.Tacq is None:
        Tacq = 5
    else:
        Tacq = args.Tacq

    print("Run acquisition : "+str(Tacq))
    run(Tacq)
