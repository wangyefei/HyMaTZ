# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 15:10:00 2017

@author: Fei
"""

import sys
try:
    from PyQt4.QtCore import (PYQT_VERSION_STR, QFile, QFileInfo, QSettings,
            QT_VERSION_STR, QTimer, QVariant, Qt, SIGNAL)
    from PyQt4.QtGui import (QAction, QActionGroup, QApplication,
            QDockWidget, QFileDialog, QFrame, QIcon, QImage, QImageReader,
            QImageWriter, QInputDialog, QKeySequence, QLabel, QListWidget,
            QMainWindow, QMessageBox, QPainter, QPixmap, QPrintDialog,QSplitter,QLineEdit,QComboBox,QCheckBox,
            QPrinter, QSpinBox,QStandardItemModel,QWidget,QGridLayout,QTableView,QDialog,QStandardItem,QPushButton)
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar    
except:
    from PyQt5.QtWidgets import (QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar     

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure



class iron_profile(QDialog):
    
    def __init__(self,Depth = None,Pressure = None, Temperature = None,usercontrol=None,fo=None,mgwa=None,mgri=None):
        super(iron_profile,self).__init__()
        self.Depth = np.array(Depth)
        self.Pressure = np.array(Pressure)
        self.Temperature = np.array(Temperature)
        self.OL_iron = fo
        self.WA_iron = mgwa
        self.RI_iron = mgri
        self.usercontrol=usercontrol
        self.OL_function = self.usercontrol['OL_function']
        self.WA_function = self.usercontrol['WA_function']
        self.RI_function = self.usercontrol['RI_function']
        
        self.figure = Figure(dpi=200)
        self.ax = self.figure.add_subplot(111)  
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
      
        
        
        layout_main = QGridLayout()
        layout_main.addWidget(self.canvas,0,0,1,1) 
        layout_main.addWidget(self.toolbar,1,0,1,1) 
        self.BTN_layout()
        self.iron_funtion()
        layout_main.addLayout(self.layout_BTN,2,0,1,1) 
        
        self.LAYOUT = QGridLayout()
        self.LAYOUT.addLayout(self.layout_control,0,0,1,1)
        self.LAYOUT.addLayout(layout_main,0,1,1,1)
        self.setLayout(self.LAYOUT)

        for i in range(len(self.Depth)):
            if self.Depth[i] >=410:
                n1=i
                break
        for i in range(n1,len(self.Depth)):
            if self.Depth[i] >=520:
                n2=i
                break 
        #print (n1,n2)        
        self.ax.plot(self.Depth[0:n1],self.OL_iron[0:n1]*100,'b')
        self.ax.plot(self.Depth[n1:n2],self.WA_iron[n1:n2]*100,'r')
        self.ax.plot(self.Depth[n2:],self.RI_iron[n2:]*100,'y')
        self.ax.set_ylim(-1,101)
        self.ax.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
        self.ax.set_ylabel('iron content (fo number)', fontname="Times New Roman",fontsize=10)
        self.canvas.draw()  
        
        
    def BTN_layout(self):
        self.layout_BTN = QGridLayout()
        

    
    def iron_funtion(self):
        
        self.layout_control =   QGridLayout()
        Label= QLabel()
        Label.setText(" Profile: change function to P,T or D \n Olivine iron function: \n Fo number = x*P + y*T + z*D")
        Label.setWordWrap(True) 
        self.layout_control.addWidget(Label,1,0,1,1)
        
        
        
        
        Label= QLabel()
        Label.setText('Olivine iron function') 
        self.layout_control.addWidget(Label,2,0,1,1)
        Label= QLabel()
        Label.setText('Wadsleyite iron function') 
        self.layout_control.addWidget(Label,4,0,1,1)
        Label= QLabel()
        Label.setText('Ringwoodite iron function') 
        self.layout_control.addWidget(Label,6,0,1,1)        
        
        self.textEdit_OL1 = QLineEdit()
        self.textEdit_OL1.setFixedWidth(200)
        self.textEdit_OL1.setText(self.OL_function)
        self.textEdit_OL1.returnPressed.connect(self.editor_OL)
        self.layout_control.addWidget(self.textEdit_OL1,3,0,1,1)
        #self.connect(self.textEdit_OL1,QtCore.SIGNAL("returnPressed()"),self.editor_OL)
        
        self.textEdit_WA1 = QLineEdit()
        self.textEdit_WA1.setFixedWidth(200)
        self.textEdit_WA1.setText(self.WA_function)
        self.textEdit_WA1.returnPressed.connect(self.editor_WA)
        self.layout_control.addWidget(self.textEdit_WA1,5,0,1,1)

        
        self.textEdit_RI1 = QLineEdit()
        self.textEdit_RI1.setFixedWidth(200)
        self.textEdit_RI1.setText(self.RI_function)
        self.textEdit_RI1.returnPressed.connect(self.editor_RI)
        self.layout_control.addWidget(self.textEdit_RI1,7,0,1,1)

        

    
        self.Button_change = QPushButton("Update")
        self.Button_change.clicked.connect(self.plot)
        self.layout_control.addWidget(self.Button_change,8,0,1,1)
        self.Button_change.setAutoDefault(False)
        
                
    def editor_OL(self):
        self.OL_function = self.textEdit_OL1.text()
        D = self.Depth;P=self.Pressure;T = self.Temperature
        self.OL_iron = (eval(self.OL_function) +np.zeros(len(D)))/100
        self.usercontrol['OL_function']=self.OL_function


        
        
    def editor_WA(self):
        self.WA_function = self.textEdit_WA1.text()
        D = self.Depth;P=self.Pressure;T = self.Temperature
        self.WA_iron = (eval(self.WA_function) +np.zeros(len(D)))/100
        self.usercontrol['WA_function']=self.WA_function

        
    def editor_RI(self):
        self.RI_function = self.textEdit_RI1.text()
        D = self.Depth;P=self.Pressure;T = self.Temperature
        self.RI_iron = (eval(self.RI_function) +np.zeros(len(D)))/100
        self.usercontrol['RI_function']=self.RI_function

   
    
    
    
    def plot(self):
        self.editor_OL()
        self.editor_WA()
        self.editor_RI()
        self.ax.clear()
        for i in range(len(self.Depth)):
            if self.Depth[i] >=410:
                n1=i
                break
        for i in range(n1,len(self.Depth)):
            if self.Depth[i] >=520:
                n2=i
                break 
        #print (n1,n2)   
        self.ax.plot(self.Depth[0:n1],self.OL_iron[0:n1]*100,'b')
        self.ax.plot(self.Depth[n1:n2],self.WA_iron[n1:n2]*100,'r')
        self.ax.plot(self.Depth[n2:],self.RI_iron[n2:]*100,'y')
        self.ax.set_ylim(0,100)
        self.ax.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
        self.ax.set_ylabel('iron content (fo number)', fontname="Times New Roman",fontsize=10)
        self.canvas.draw()  
    
    def RETURN(self):
        return [self.OL_iron,self.WA_iron,self.RI_iron]
    

        
if __name__ == "__main__": 
    n=300
    D=np.linspace(0,800,300);P=np.zeros(n);T=np.zeros(n);
    usercontrol = {
            'OL_function':'4.0*1e-6*D*D',
            'WA_function':'3.3',
            'RI_function':'1.7',
            }
    qApp = QApplication(sys.argv)
    aw = iron_profile(D,P,T,usercontrol,np.ones(300),np.ones(300),np.ones(300))
    aw.show()
    sys.exit(qApp.exec_())
    print (aw.usercontrol['OL_function'])