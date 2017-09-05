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



class water_profile(QDialog):
    
    def __init__(self,Depth = None,Pressure = None, Temperature = None,usercontrol=None):
        super(water_profile,self).__init__()
        self.Depth = np.array(Depth)
        self.Pressure = np.array(Pressure)
        self.Temperature = np.array(Temperature)
        self.OL_water = np.ones(len(Depth))
        self.WA_water = np.zeros(len(Depth))
        self.RI_water = np.zeros(len(Depth))
        self.usercontrol=usercontrol
        self.OL_function = self.usercontrol['OL_function']
        self.WA_function = self.usercontrol['WA_function']
        self.RI_function = self.usercontrol['RI_function']
        self.OL=self.usercontrol['OLcoefficient']
        self.WA=self.usercontrol['WAcoefficient']
        self.RI=self.usercontrol['RIcoefficient']
        
        self.figure = Figure(dpi=200)
        self.ax = self.figure.add_subplot(111)  
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
      
        
        
        layout_main = QGridLayout()
        layout_main.addWidget(self.canvas,0,0,1,1) 
        layout_main.addWidget(self.toolbar,1,0,1,1) 
        self.BTN_layout()
        self.water_funtion()
        layout_main.addLayout(self.layout_BTN,2,0,1,1) 
        
        self.LAYOUT = QGridLayout()
        self.LAYOUT.addLayout(self.layout_control,0,0,1,1)
        self.LAYOUT.addLayout(layout_main,0,1,1,1)
        self.setLayout(self.LAYOUT)
        self.plot()
        
        
    def BTN_layout(self):
        self.layout_BTN = QGridLayout()

        

        Label= QLabel()
        Label.setText('Olivine water content (% Max capacity)') 
        self.layout_BTN.addWidget(Label,0,0,1,1)
        Label= QLabel()
        Label.setText('Wadsleyite water content (% Max capacity)') 
        self.layout_BTN.addWidget(Label,0,1,1,1)
        Label= QLabel()
        Label.setText('Ringwoodite water content (% Max capacity)') 
        self.layout_BTN.addWidget(Label,0,2,1,1)
        
        self.textEdit_OL =  QSpinBox()
        self.textEdit_OL.setMaximum(100)
        self.textEdit_OL.setMinimum(0)
        self.textEdit_OL.setValue(self.OL)
        self.textEdit_OL.valueChanged.connect(self.editor_OL)        
        self.layout_BTN.addWidget(self.textEdit_OL,1,0,1,1)

        self.textEdit_WA = QSpinBox()
        self.textEdit_WA.setMaximum(100)
        self.textEdit_WA.setMinimum(0)
        self.textEdit_WA.setValue(self.WA)
        self.textEdit_WA.valueChanged.connect(self.editor_WA)          
        self.layout_BTN.addWidget(self.textEdit_WA,1,1,1,1)
       
        self.textEdit_RI = QSpinBox()
        self.textEdit_RI.setMaximum(100)
        self.textEdit_RI.setMinimum(0)
        self.textEdit_RI.setValue(self.RI)
        self.textEdit_RI.valueChanged.connect(self.editor_RI)          
        self.layout_BTN.addWidget(self.textEdit_RI,1,2,1,1)
        
    
    def water_funtion(self):
        
        self.layout_control =   QGridLayout()
        Label= QLabel()
        Label.setText(" Profile: change function to P,T or D \n Olivine water function: \n H₂O (wt%) = x*P + y*T + z*D")
        Label.setWordWrap(True) 
        self.layout_control.addWidget(Label,1,0,1,1)
        
        
        
        
        Label= QLabel()
        Label.setText('Olivine water function') 
        self.layout_control.addWidget(Label,2,0,1,1)
        Label= QLabel()
        Label.setText('Wadsleyite water function') 
        self.layout_control.addWidget(Label,4,0,1,1)
        Label= QLabel()
        Label.setText('Ringwoodite water function') 
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
        self.OL=float(self.textEdit_OL.value()) 
        self.OL_function = self.textEdit_OL1.text()
        D = self.Depth;P=self.Pressure;T = self.Temperature
        self.OL_water = (eval(self.OL_function) +np.zeros(len(D)))*self.OL/100
        for i in range(len(D)):
            if self.OL_water[i] >=1.1:
                self.OL_water[i] = 1.1
        self.usercontrol['OL_function']=self.OL_function
        self.usercontrol['OLcoefficient'] = self.OL

        
        
    def editor_WA(self):
        self.WA=float(self.textEdit_WA.value())
        self.WA_function = self.textEdit_WA1.text()
        D = self.Depth;P=self.Pressure;T = self.Temperature
        self.WA_water = (eval(self.WA_function) +np.zeros(len(D)))*self.WA/100
        for i in range(len(D)):
            if self.WA_water[i] >=3.3:
                self.WA_water[i] = 3.3
        self.usercontrol['WA_function']=self.WA_function
        self.usercontrol['WAcoefficient'] = self.WA
        
    def editor_RI(self):
        self.RI=float(self.textEdit_RI.value()) 
        self.RI_function = self.textEdit_RI1.text()
        D = self.Depth;P=self.Pressure;T = self.Temperature
        self.RI_water = (eval(self.RI_function) +np.zeros(len(D)))*self.RI/100
        for i in range(len(D)):
            if self.RI_water[i] >=2.7:
                self.RI_water[i] = 2.7
        self.usercontrol['RI_function']=self.RI_function
        self.usercontrol['RIcoefficient'] = self.RI
   
    
    
    
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
        self.ax.plot(self.Depth[0:n1],self.OL_water[0:n1],'b')
        self.ax.plot(self.Depth[n1:n2],self.WA_water[n1:n2],'r')
        self.ax.plot(self.Depth[n2:],self.RI_water[n2:],'y')
        self.ax.set_ylim(0,3.5)
        self.ax.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
        self.ax.set_ylabel('Water content (wt%)', fontname="Times New Roman",fontsize=10)
        self.canvas.draw()  
    
    def RETURN(self):
        return [self.OL_water,self.WA_water,self.RI_water]
    

        
if __name__ == "__main__": 
    n=300
    D=np.linspace(0,800,300);P=np.zeros(n);T=np.zeros(n);
    usercontrol = {
            'OL_function':'4.0*1e-6*D*D',
            'WA_function':'3.3',
            'RI_function':'1.7',
            'OLcoefficient':10,
            'WAcoefficient':20,
            'RIcoefficient':30,
            }
    qApp = QApplication(sys.argv)
    aw = water_profile(D,P,T,usercontrol)
    aw.show()
    sys.exit(qApp.exec_())
    print (aw.usercontrol['OL_function'])