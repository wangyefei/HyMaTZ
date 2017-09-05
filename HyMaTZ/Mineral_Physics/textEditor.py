'''
This class provide a GUI text editor. 
'''
import sys
import os
try:
    from PyQt4.QtGui import (QWidget, QPushButton, QLineEdit, QLabel,QDialog,
                                 QGridLayout, QApplication,QSpinBox,QTextEdit)
except:
    from PyQt5.QtWidgets import (QWidget, QPushButton,QDialog,
                                 QGridLayout, QApplication,QTextEdit)


class TextEditor(QDialog):
    '''
    This is class is a basic class use to show txt file and save data
    '''
    
    def __init__(self,name=None):
        super(TextEditor,self).__init__()
        self.layout = QGridLayout()
        self.name = name
        self.text = QTextEdit(self)
        self.set_Text()
        self.layout.addWidget(self.text)
        self.save_BTN = QPushButton('save')
        self.save_BTN.clicked.connect(self.save)
        self.layout.addWidget(self.save_BTN)        
        self.setLayout(self.layout)     
        self.resize(500,500)
    
    def save(self):
        f = open(self.address, 'w')
        filedata = self.text.toPlainText()
        f.write(filedata)   
        f.close()
            
            
    def set_Text(self):
        try:
            address=os.path.join(os.path.dirname(__file__),'EXPDATA',self.name+'.txt')
            self.address=address
        except:
            address=os.path.join(os.path.dirname(__file__),'Mineral','EXPDATA',self.name+'.txt')
            self.address=address
            
        f = open(self.address, 'r')
        filedata = f.read()
        self.text.setText(filedata)     
        f.close()

        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    main= TextEditor(name='Olivine')
    main.show()
    sys.exit(app.exec_())