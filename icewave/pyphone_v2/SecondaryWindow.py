
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



class SecondaryWindow(QWidget):
    def __init__(self):
        super().__init__()        
        self.pagelayout = QVBoxLayout()


    def add_buttons(self,buttons):
        self.button_layout = QHBoxLayout()
        self.pagelayout.addLayout(self.button_layout)
      
        for key in buttons.keys():
            btn = QPushButton(key)
            btn.pressed.connect(buttons[key])
            self.button_layout.addWidget(btn)
        
    def add_console(self):
        self.console = Color('black')
        self.console.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.console.text=""
        self.console.setText(self.console.text)           
        self.pagelayout.addWidget(self.console)
    
    def pt(self,s,keep=True):
        print("it is printing!")
        s = str(s)
        if keep:
            s=self.console.text+"\n"+s
        self.console.setText(s)
        self.mainWindow.repaint()
        self.console.text=s
        print('over')

 
    def base_window(self,buttons={},color='black',text='ok'):
        pagelayout = QVBoxLayout()
        button_layout = QHBoxLayout()
        pagelayout.addLayout(button_layout)
        
        for key in buttons.keys():
            btn = QPushButton(key)
            btn.pressed.connect(buttons[key])
            button_layout.addWidget(btn)
        
        pagelayout.addWidget(Color('black'))

        widget = QWidget()
        widget.setLayout(pagelayout)
        #self.console.label.setText("Console is here")
        return widget

    def parameter_box(self,label,default=10):
        wid1 = QWidget()
        
        layout = QHBoxLayout()
        layout.addWidget(QLabel(label))
        
        wid2 = QLineEdit()
        wid2.setMaxLength(10)
        wid2.setText(str(default))
        #widget.setReadOnly(True) # uncomment this to make readonly

#        wid2.returnPressed.connect(self.return_pressed)
#        wid2.selectionChanged.connect(self.selection_changed)
        #wid.textChanged.connect(lambda x:self.text_changed(x,var))
#        wid2.textEdited.connect(self.text_edited)   
        layout.addWidget(wid2)
        wid1.setLayout(layout)
        return wid1,wid2


#    def add_parameter_box(self,labels,defaults=[30]):
#   	    parameter_layout = QHBoxLayout()
        #self.pagelayout.addLay
#        widget.pagelayout.addLayout(parameter_layout)

"""    
    def add_parameter_box(self,labels,defaults=[30]):
   	    parameter_layout = QHBoxLayout()

        for label in labels:
           label = 'Acquisition time (s) :'
            self.temperature_param_Tacq = 30

            wid1,wid2 = self.parameter_box(label)
            wid2.textChanged.connect(self.temperature_Tacq)

            parameter_layout.addWidget(wid1)
"""