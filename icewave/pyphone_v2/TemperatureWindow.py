

import sys

import subprocess
import time
import numpy as np

from pprint import pprint

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
import battery_T
import connect
#import run_adb_usb as adb_usb 
#import game_breaker  
#import run_phyphox 


class TemperatureWindow(SecondaryWindow):
    def __init__(self,mainWindow):
        super().__init__()        
#        self.pagelayout = QVBoxLayout()
        self.mainWindow = mainWindow
        buttons = {'Run':self.act_temperature_run,'Display':self.act_temperature_display}        
        self.add_buttons(buttons)
        print(self.pagelayout.count())
        
        parameter_layout = QHBoxLayout()
        self.pagelayout.addLayout(parameter_layout)

        label = 'Acquisition time (s) :'
        self.param_Tacq = 10

        wid1,wid2 = self.parameter_box(label,default=self.param_Tacq)
        wid2.textChanged.connect(self.temperature_Tacq)
        
        parameter_layout.addWidget(wid1)
        
        self.gridlayout = QGridLayout()
        self.pagelayout.addLayout(self.gridlayout)
                
        adrlist,self.phonelist = connect.connect() 
        
        self.phoneboxes={}
        for phone in self.phonelist:
            i = int(np.floor(phone/10))
            j = phone-i*10
            print(i,j)
            
            self.gridlayout.addWidget(self.click_buttons(phone),i,j)
            #self.gridlayout.addChildWidget(widget)
            
        print('number of checkbox:'+str(self.gridlayout.count()))
        
        self.add_console()
        self.setLayout(self.pagelayout)
        #self.console.label.setText("Console is here")
        
    def act_temperature_run(self):
        self.pt("Run temperature measurement",keep=False)

        phone_to_run=[]
        for i,phone in enumerate(self.phonelist):
            checked = (self.phoneboxes[phone].checkState() == Qt.CheckState.Checked)
            if checked:
                print(phone)
                phone_to_run.append(phone)
        self.pt(phone_to_run,keep=False)
#            print(phone,)
#        adrlist,phonelist = battery_T.connect()
        phone_dict = battery_T.make_phone_dict(phone_to_run)

        t1 = Process(target=self.run_temperature, args=(phone_dict,))
        t1.start()

#        t1 = Process(target=battery_T.acquire, args=(phone_dict,),kwargs=keywords)
    
        #battery_T.acquire(phone_dict,Tmax=100)
            #checked = Qt.CheckState(self.gridlayout.itemAt(i)) == Qt.CheckState.Checked
            #print(phone,checked)

    def run_temperature(self,phone_dict):
        print(self.param_Tacq)
        keywords = {'Tmax': self.param_Tacq}#,'window':self}
        battery_T.acquire(phone_dict,**keywords)
        self.pt('Done')
        

    def act_temperature_display(self):
        self.pt('Display temperature')
        
    def click_buttons(self,phone):
        self.phoneboxes[phone] = QCheckBox(str(phone))
        self.phoneboxes[phone].setCheckState(Qt.CheckState.Checked)

        # For tristate: widget.setCheckState(Qt.PartiallyChecked)
        # Or: widget.setTristate(True)
        self.phoneboxes[phone].stateChanged.connect(self.mainWindow.show_state)
        return self.phoneboxes[phone]
        

    def temperature_Tacq(self, s):
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
