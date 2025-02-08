import sys

import subprocess
import time
import numpy as np

from pprint import pprint

from multiprocessing import Process
    
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
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

import connect
        
from SecondaryWindow import SecondaryWindow
from GobannosWindow import GobannosWindow

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Phone-Fleet ")

        self.tabs = QTabWidget()
        self.tabs.setTabPosition(QTabWidget.TabPosition.West)
        self.tabs.setMovable(True)

        tablist = self.get_tabs()
        for key in tablist.keys():
            #widget = self.base_window(buttons=tablist[key],text='toto')
            widget = SecondaryWindow()
            widget.add_buttons(tablist[key])
            widget.add_console()
            widget.setLayout(widget.pagelayout)
            self.tabs.addTab(widget, key)
#        tabs.addTab(widget_network,'Network')

        #self.TemperatureWindow = TemperatureWindow(self)
        #self.tabs.addTab(self.TemperatureWindow, 'Temperature')

        self.GobannosWindow=GobannosWindow(self)
        self.tabs.addTab(self.GobannosWindow,'Phyphox')
        self.setCentralWidget(self.tabs)
        
    def get_tabs(self):
        tablist = { 'Network':  {'Display':self.act_test_connect,\
                                'Reset wifi':self.act_wifi_connect}, \
                    'Time':     {'sync':self.act_time_sync,\
                                'display':self.act_time_display}, \
                    'Files':{'Explorer':self.file_explorer},\
                    'Location':{'Location':self.act_get_location}}
#                   'adb Screen':,\
#                   'Plots':[]}
        return tablist
        
    def show_state(self, s):
        print(Qt.CheckState(s) == Qt.CheckState.Checked)

    def pt(self,s):
        s = str(s)
        self.tabs.currentWidget().console.setText(s)
        self.repaint()

    def act_get_location(self):
        print("run location on Phyphox")

    def act_adb_connect(self):
        self.pt('running ...')       
        adb_usb.run()
        self.pt('done')

#        self.stacklayout.setCurrentIndex(0)

    def act_test_connect(self):
        c=subprocess.run(['adb','devices'],text=True,capture_output=True)
        log = c.stdout#.split('\n')
        #self.pt('toto')
        self.pt(log)
       
    def act_wifi_connect(self):
        base = '/home/turbots/Documents/Codes/Python/pyphone/'
        subprocess.call(base+'startup_wifi.sh')
        self.pt('Reset wifi')

    def act_time_sync(self):        
        adrlist,phonelist = connect.connect()        
        multi_line = str(phonelist) +'\n'+"Number of phone connected :"+str(len(phonelist))+'\nRunning ...'
        self.pt(multi_line)
        adb_usb.run()
        self.pt('done')
        self.table = game_breaker.compute_timetable(adrlist,phonelist)   
        #self.pt(self.table)
                
    def act_time_display(self):
        if hasattr(self,"table"):
            self.pt(self.table)
#            for tab in self.table:
#                print(int(tab[0]),float(tab[1].split('\n')[0]))
        else:
            self.pt('time table not yet defined')

    def act_temperature_display(self):
        self.pt('Measure temperature')
        
    def file_explorer(self):
        self.pt('Start file_explorer')

    def act_take_picture(self):
        phone = 0
        self.pt('take_picture with phone '+str(phone))

        t1 = Process(target=self.take_picture,args=(phone,))
        t1.start()
        
    def take_picture(self,phone):
        connect_phone.launch_camera(phone)
        time.sleep(1)
        connect_phone.take_picture(phone)


if __name__=='__main__':
                
    app = QApplication(sys.argv)

    w = MainWindow()
    w.resize(800,600);
    w.show()

    app.exec()
