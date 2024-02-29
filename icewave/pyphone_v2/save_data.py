

import time
from multiprocessing import Process

import icewave.pyphone_v2.connect as connect
import icewave.pyphone_v2.run_phyphox as run_phyphox


def save():
    adrlist,phonelist = connect.connect()  
    base = connect.ipbase()
    iplist = [base + str(100+phone) for phone in phonelist]

    print('Save data Phyphox on '+str(phonelist))
    tstamp = int(time.time())
    keywords = {'folder':'PhyPhox_'+str(tstamp)}
    
    print(iplist[0])
    url = run_phyphox.get_base_url(iplist[0])
       
       
    date = '121223_'+str(time.time())
    out = run_phyphox.save_single(url,folder=date)
    
    return out
    #battery_T.acquire(phone_dict,Tmax=100)
        #checked = Qt.CheckState(self.gridlayout.itemAt(i)) == Qt.CheckState.Checked
        #print(phone,checked)

if __name__ =='__main__':
    save()
