# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:11:37 2017

@author: Fei
"""

try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar  
    PYQT = 4
except:
    from PyQt5.QtWidgets import (QProgressDialog,QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QHBoxLayout,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (QThread,pyqtSignal,PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar     
    PYQT = 5
    
from Mineral_Physics.attenuation import Anelastic


class Attenuation_Model():
    
    def __init__(self,Name='GoesQ1',A=None,w=None,a=None,H=None,V=None,g1=None,g2=None,Type=1):
        self.Name = Name
        self.Type=Type
        self.A = A
        self.w = w
        self.a = a
        self.H = H
        self.V = V
        self.g1 = g1
        self.g2 = g2
 
GoesQ1 = Attenuation_Model(Name = 'Goes et al.(2000) Q1',A=1.48*1e-1,w=0.05,a=0.15,H=500,V=20)    
GoesQ2 = Attenuation_Model(Name = 'Goes et al.(2000) Q2',A=2.0*1e-4,w=0.05,a=0.25,H=584,V=21)  
CammaranoQ1 = Attenuation_Model(Name = 'Cammarano et al.(2003) Q1',A=0.5,w=1,a=0.2,g1=20,g2=10,Type=0)
CammaranoQ2 = Attenuation_Model(Name = 'Cammarano et al.(2003) Q2',A=0.8,w=1,a=0.2,g1=20,g2=10,Type=0)
CammaranoQ3 = Attenuation_Model(Name = 'Cammarano et al.(2003) Q3',A=1.1,w=1,a=0.2,g1=20,g2=10,Type=0)
CammaranoQ4 = Attenuation_Model(Name = 'Cammarano et al.(2003) Q4',A=0.035,w=1,a=0.2,g1=30,g2=15,Type=0)
CammaranoQ5 = Attenuation_Model(Name = 'Cammarano et al.(2003) Q5',A=0.056,w=1,a=0.2,g1=30,g2=15,Type=0)
model_list = [GoesQ1,GoesQ2,CammaranoQ1,CammaranoQ2,CammaranoQ3,CammaranoQ4,CammaranoQ5]
    
    
class Attenuation(QDialog):
    
    def __init__(self,model_list=model_list,parameter_list=[2.0*1e-4,0.05,0.25,584,21,0,0]):
        super(Attenuation,self).__init__()
        layout = QGridLayout()
        self.model_list=model_list
        
        i=0;j=0
        self.BTN_Attenuation=QComboBox()
        self.BTN_Attenuation.addItem("Select a Model")
        for k in range(len(self.model_list)):
            self.BTN_Attenuation.addItem(self.model_list[k].Name)
        layout.addWidget(self.BTN_Attenuation,i,0,1,6)
        self.BTN_Attenuation.currentIndexChanged.connect(self.Model_select)        
        
        i+=1
        self.Label1 = QLabel("A:");
        self.textEdit1 = QLineEdit()
        self.textEdit1.setText(str(parameter_list[0]))
        self.textEdit1.returnPressed.connect(self.editor1)
        layout.addWidget(self.Label1,i,j,1,1);layout.addWidget(self.textEdit1,i+1,j,1,1)
        self.Label2 = QLabel("w (Hz):");j+=1;
        self.textEdit2 = QLineEdit()
        self.textEdit2.setText(str(parameter_list[1]))
        self.textEdit2.returnPressed.connect(self.editor2)
        layout.addWidget(self.Label2,i,j,1,1) ;layout.addWidget(self.textEdit2,i+1,j,1,1)         
        self.Label3 = QLabel("a:");j+=1;
        self.textEdit3 = QLineEdit()
        self.textEdit3.setText(str(parameter_list[2]))
        self.textEdit3.returnPressed.connect(self.editor3)
        layout.addWidget(self.Label3,i,j,1,1) ;layout.addWidget(self.textEdit3,i+1,j,1,1)             
        self.Label4 = QLabel("H (KJ/mol):");j+=1;
        self.textEdit4 = QLineEdit()
        self.textEdit4.setText(str(parameter_list[3]))
        self.textEdit4.returnPressed.connect(self.editor4)
        layout.addWidget(self.Label4,i,j,1,1)  ;layout.addWidget(self.textEdit4,i+1,j,1,1)  
        self.Label5 = QLabel("V (cm3/mol):");j+=1;
        self.textEdit5 = QLineEdit()
        self.textEdit5.setText(str(parameter_list[4]))
        self.textEdit5.returnPressed.connect(self.editor5)
        layout.addWidget(self.Label5,i,j,1,1);layout.addWidget(self.textEdit5,i+1,j,1,1)
        self.Label6 = QLabel("g1:");j+=1;
        self.textEdit6 = QLineEdit()
        self.textEdit6.setText(str(parameter_list[5]))
        self.textEdit6.returnPressed.connect(self.editor6)
        layout.addWidget(self.Label6,i,j,1,1);layout.addWidget(self.textEdit6,i+1,j,1,1)
        self.Label7 = QLabel("g2:");j+=1;
        self.textEdit7 = QLineEdit()
        self.textEdit7.setText(str(parameter_list[6]))
        self.textEdit7.returnPressed.connect(self.editor7)
        layout.addWidget(self.Label7,i,j,1,1);layout.addWidget(self.textEdit7,i+1,j,1,1)
        
        i+=2
        self.okButton =QPushButton("&OK",self)
        self.okButton.clicked.connect(self.OK)
        self.okButton.setAutoDefault(False)
        #self.cancelButton = QPushButton("Cancel",self)
        #self.cancelButton.clicked.connect(self.Cancel)
        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.okButton)
        #buttonLayout.addWidget(self.cancelButton)          
        layout.addLayout(buttonLayout, i, 5, 1, 2)
        
        
        self.setLayout(layout)
        self.setWindowTitle("Attenuation")
        self.resize(100,50)        


    def editor1(self):
        self.A = (float(self.textEdit1.text()))
    def editor2(self):
        self.w = (float(self.textEdit2.text()))
    def editor3(self):
        self.a = (float(self.textEdit3.text()))        
    def editor4(self):
        self.H = (float(self.textEdit4.text()))
    def editor5(self):
        self.V = (float(self.textEdit5.text()))
    def editor6(self):
        self.g1 = (float(self.textEdit6.text()))
    def editor7(self):
        self.g2 = (float(self.textEdit7.text()))
        
    def  Model_select(self,a):
        if a == 0:
            pass
        else:
            model = self.model_list[a-1]
            if model.Type==1:
                self.Type=1
                self.textEdit4.setDisabled(False)
                self.textEdit5.setDisabled(False)
                self.textEdit1.setText(str(model.A))
                self.textEdit2.setText(str(model.w))
                self.textEdit3.setText(str(model.a))
                self.textEdit4.setText(str(model.H))
                self.textEdit5.setText(str(model.V))
                self.textEdit6.setText('')
                self.textEdit6.setDisabled(True)
                self.textEdit7.setDisabled(True)
            if model.Type==0:
                self.Type=0
                self.textEdit6.setDisabled(False)
                self.textEdit7.setDisabled(False)                
                self.textEdit1.setText(str(model.A))
                self.textEdit2.setText(str(model.w))
                self.textEdit3.setText(str(model.a))
                self.textEdit4.setText('');self.textEdit4.setDisabled(True)
                self.textEdit5.setText('');self.textEdit5.setDisabled(True)
                self.textEdit6.setText(str(model.g1))
                self.textEdit7.setText(str(model.g2))
       

       
    def OK(self):
        self.editor1()
        self.editor2()
        self.editor3()        
        if self.Type==1:
            self.editor4()
            self.editor5()
        if self.Type==0:
            self.editor6()  
            self.editor7()  
        self.close()
    def Cancel(self):
        self.close()        


if __name__ == "__main__":
    import sys
    qApp = QApplication(sys.argv)

    aw = Attenuation()
    aw.show()
    sys.exit(qApp.exec_())  