
import sys

import subprocess
import time
import numpy as np

from pprint import pprint,pformat

from multiprocessing import Process

    
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QStackedLayout,
    QVBoxLayout,
    QGridLayout,
    QTabWidget,
    QWidget,
    QCheckBox,
    QLineEdit,
    QTextEdit
)

from layout_colorwidget import Color

from SecondaryWindow import SecondaryWindow
import connect
import run_phyphox

class PhyphoxWindow(SecondaryWindow):
    def __init__(self,mainWindow):
        super().__init__()        
        self.pagelayout = QVBoxLayout()
        self.mainWindow = mainWindow


        buttons_list = [{},{}]
        buttons_list[0] = { 'Start Phyphox':self.act_phyphox_start,\
                            'Check':self.act_phyphox_check,\
                            'Config':self.act_phyphox_config,\
                            'Location':self.act_phyphox_location}                        
        buttons_list[1] = {'Run':self.act_phyphox_run,\
                    'Clear':self.act_phyphox_clear,\
                    'Stop':self.act_phyphox_stop,\
                    'Save':self.act_phyphox_save}
        
        buttons_layout = QVBoxLayout()
        self.pagelayout.addLayout(buttons_layout)
        
        for buttons in buttons_list:  
            button_layout = QHBoxLayout()
            buttons_layout.addLayout(button_layout)
            for key in buttons.keys():
                btn = QPushButton(key)
                btn.pressed.connect(buttons[key])
                button_layout.addWidget(btn)
        
        parameter_layout = QHBoxLayout()
        buttons_layout.addLayout(parameter_layout)

        label = 'Acquisition time (s) :'
        self.param_Tacq = 5

        wid1,wid2 = self.parameter_box(label,default=self.param_Tacq)
        wid2.textChanged.connect(self.set_Tacq)
        
        parameter_layout.addWidget(wid1)
        
        self.layout_Phy = QGridLayout()
        self.pagelayout.addLayout(self.layout_Phy)
        
        adrlist,self.phonelist = connect.connect()   
        
        self.phoneboxes={}
        for phone in self.phonelist:
            self.add_box(phone)
            
        print('number of checkbox:'+str(self.layout_Phy.count()))

        self.add_console()
        self.setLayout(self.pagelayout)
        #self.console.label.setText("Console is here")

    def add_box(self,phone):
        i = int(np.floor(phone/10))
        j = phone-i*10
        print(i,j)    
        self.layout_Phy.addWidget(self.click_button_Phy(phone),i,j)
        #self.gridlayout.addChildWidget(widget)
                
    def add_buttons(self,buttons):
        self.button_layout = QHBoxLayout()
        self.pagelayout.addLayout(self.button_layout)
      
        for key in buttons.keys():
            btn = QPushButton(key)
            btn.pressed.connect(buttons[key])
            self.button_layout.addWidget(btn)
        
    def set_Tacq(self, s):
        print("Text changed...")
        if len(s)>0:
            try:
                self.param_Tacq=int(s)
            except:
                self.pt("not an integer number")
                self.param_Tacq=0
        else:
            self.param_Tacq=0
        print(self.param_Tacq)

    def click_button_Phy(self,phone):
        self.phoneboxes[phone] = QCheckBox(str(phone))
        self.phoneboxes[phone].setCheckState(Qt.CheckState.Checked)#unchecked
        self.phoneboxes[phone].stateChanged.connect(self.mainWindow.show_state)
        return self.phoneboxes[phone]        
        
    def act_phyphox_config(self):
        adrlist,phonelist = connect.connect()  
        base = connect.ipbase()
        iplist = [base + str(100+phone) for phone in phonelist]
   
        self.pt(iplist,keep=False)
        R = run_phyphox.run_config(iplist)
        s = pformat(R[0]['buffers'])
        self.pt(s)
        #t1 = Process(target=run_phyphox.run_config, args=(iplist,))
        #t1.start()
   
    def act_phyphox_location(self):
        adrlist,phonelist = connect.connect()  
        base = connect.ipbase()
        iplist = [base + str(100+phone) for phone in phonelist]
        self.pt("Get location from Phyphox",keep=False)
        locations = run_phyphox.run_get(iplist,'location')
        
        self.display_locations(locations)        
        #t1 = Process(target=run_phyphox.run_config, args=(iplist,))
        #t1.start()

    def display_locations(self,locations):
        s = ""
        adrlist,phonelist = connect.connect()  

        print(locations)
        for phone in phonelist:
            r = locations[phone]
            s = s+str(phone)+'\t'
            s = s+'Longitude : '+str(r['buffer']['locLon']['buffer'][0])+'\t'
            s = s+'Latitude : '+str(r['buffer']['locLat']['buffer'][0])
            s = s+'\n'
        print(s)
        self.pt(s)


    def act_phyphox_start(self):
        self.act_phyphox_check(self)

        self.pt('starting phyphox on '+str(self.phone_toconnect))
        for phone in self.phone_toconnect:
            connect.launch_phyphox(phone)            
        self.act_phyphox_check(self)
        
    def act_phyphox_clear(self):
        adrlist,self.phonelist = connect.connect()   
        self.pt('clear phyphox',keep=False)
    
    def act_phyphox_check(self):
        adrlist,self.phonelist = connect.connect()   
        self.pt("Checking Phyphox phone connexions, please wait ...",keep=False)

        self.phyphox_check()       
#        t1 = Process(target=self.phyphox_check)
#        t1.start()
        
    def phyphox_check(self):
        self.phone_phyphox = []
        self.phone_toconnect = []

        print(self.phonelist)
        iplist = connect.get_adresslist(self.phonelist)
        R = run_phyphox.run_fun(run_phyphox.clear,iplist)
        pprint(R)
        for phone in self.phonelist:
            if R[phone]['result']:#connexion succesfull
                if phone not in self.phone_phyphox:
                    self.phone_phyphox.append(phone)
            else:
                if phone in self.phone_phyphox:
                    self.phone_phyphox.remove(phone)
                self.phone_toconnect.append(phone)

        self.pt("Phones available : "+str(self.phone_phyphox))
        self.pt("Phones to connect : "+str(self.phone_toconnect))

        self.refresh_Checkbox_Phy()
            
    def refresh_Checkbox_Phy(self):
        for phone in self.phone_phyphox:
            if not phone in self.phoneboxes.keys():
                self.add_box(phone)  
        self.show()
        
    def act_phyphox_run(self):
        phone_to_run=[]
        for i,phone in enumerate(self.phonelist):
            checked = (self.phoneboxes[phone].checkState() == Qt.CheckState.Checked)
            if checked:
                print(phone)
                phone_to_run.append(phone)
        #self.pt(phone_to_run)
        
        base = connect.ipbase()
        iplist = [base + str(100+phone) for phone in phone_to_run]
    
        self.pt('Run Phyphox on '+str(phone_to_run))
        tstamp = int(time.time())
        keywords = {'T': self.param_Tacq,'folder':'PhyPhox_'+str(tstamp)}
        
        t1 = Process(target=run_phyphox.run_serie, args=(iplist,),kwargs=keywords)
        t1.start()
        #battery_T.acquire(phone_dict,Tmax=100)
            #checked = Qt.CheckState(self.gridlayout.itemAt(i)) == Qt.CheckState.Checked
            #print(phone,checked)

    def act_phyphox_stop(self):
        adrlist,self.phonelist = connect.connect()
        self.pt('save phyphox stop')
    
    
    def act_phyphox_save(self):
        adrlist,self.phonelist = connect.connect()   
        self.pt('save phyphox data')