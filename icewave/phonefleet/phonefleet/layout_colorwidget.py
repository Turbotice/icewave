from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QPalette
from PyQt6.QtWidgets import QWidget,QLabel,QVBoxLayout,QTextEdit

class Color(QTextEdit):
    def __init__(self, color):
        super().__init__()
        #self.setAutoFillBackground(True)
        self.setReadOnly(True) # uncomment this to make readonly
        
        
        self.setTextBackgroundColor(QColor("black"))
        self.setTextColor(QColor("white"))
        self.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.setText('')   
        
        palette = self.viewport().palette()
        palette.setColor(self.viewport().backgroundRole(), QColor("black"))
#        palette.setColor(QPalette.ColorRole.Text, QColor("black"))
        self.viewport().setPalette(palette)
        

        
#        self.label = QLabel('')
 #       self.label.setPalette(palette)
 #       self.label.setAlignment(Qt.AlignmentFlag.AlignLeft)
 #       self.label.setWordWrap(True)
 #       vbox = QVBoxLayout()
 #       vbox.addWidget(self.label)
#       self.setLayout(vbox)
