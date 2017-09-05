# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 16:12:04 2017

@author: Fei
"""

import os
import numpy as np

try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *    
except:
    from PyQt5.QtWidgets import (QAbstractItemView,QVBoxLayout,QSizePolicy,QProgressDialog,QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QHBoxLayout,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (QThread,pyqtSignal,PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt) 



from Mineral_Physics.Thermodymic_model_MP import Harzburgite
from Mineral_Physics.Thermodymic_model_Phase_diagram import Harzburgite_Phase

class Layers():
    
    def __init__(self,layer_name=None,layer_depth=None,n=240000,index_P=600,index_T=400,):
        self.layer=[]
        self.layer_phase=[]
        #self.address = address
        for i in layer_name:
            self.layer.append(Harzburgite([0,0,0,0,0,0],i,n=n,index_P=index_P,index_T=index_T,address = os.path.dirname(os.path.realpath(__file__))))
            self.layer_phase.append(Harzburgite_Phase([0,0,0,0,0,0],i,n=n,index_P=index_P,index_T=index_T,address = os.path.dirname(os.path.realpath(__file__))))
        self.layer_depth = layer_depth
        
    
class GUI_layers(QDialog):
    
    def __init__(self,Depth=[],Phase=[],ComboNumber=[]):
        super(GUI_layers,self).__init__()
        self.Change=True
        self.Depth=Depth
        self.Phase=Phase
        self.CombosNumber=ComboNumber
        self.Layout()
        self.Reopen()

    def Layout(self):
        i=0
        layout = QGridLayout()
        label1 = QLabel("Layer:")
        label3 = QLabel("Depth End (km,last one must be 800):")
        label4 = QLabel("Model Type:")
        layout.addWidget(label1,i,0,1,1)   
        layout.addWidget(label3,i,1,1,1)   
        layout.addWidget(label4,i,2,1,2)   
        
        i+=1
        label1 = QLabel("Layer1:")
        self.Depth1 = QLineEdit("")
        self.BTN_Mantle_cambo1=QComboBox()
        self.BTN_Mantle_cambo1.addItem("Select a Model")
        self.BTN_Mantle_cambo1.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo1.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo1.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo1.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo1.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo1.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo1.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo1.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo1.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo1.addItem('Pyrolite5 (Ringwood, 1962*) ')  
        layout.addWidget(label1,i,0,1,1)  
        layout.addWidget(self.Depth1,i,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo1,i,2,1,2)
        
        i+=1
        label1 = QLabel("Layer2:")
        self.Depth2 = QLineEdit("")
        self.BTN_Mantle_cambo2=QComboBox()
        self.BTN_Mantle_cambo2.addItem("Select a Model")
        self.BTN_Mantle_cambo2.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo2.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo2.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo2.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo2.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo2.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo2.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo2.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo2.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo2.addItem('Pyrolite5 (Ringwood, 1962*) ')  
        layout.addWidget(label1,i,0,1,1)  
        layout.addWidget(self.Depth2,i,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo2,i,2,1,2)


        i+=1
        label1 = QLabel("Layer3:")
        self.Depth3 = QLineEdit("")
        self.BTN_Mantle_cambo3=QComboBox()
        self.BTN_Mantle_cambo3.addItem("Select a Model")
        self.BTN_Mantle_cambo3.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo3.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo3.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo3.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo3.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo3.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo3.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo3.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo3.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo3.addItem('Pyrolite5 (Ringwood, 1962*) ')  
        layout.addWidget(label1,i,0,1,1)  
        layout.addWidget(self.Depth3,i,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo3,i,2,1,2)


        i+=1
        label1 = QLabel("Layer4:")
        self.Depth4 = QLineEdit("")
        self.BTN_Mantle_cambo4=QComboBox()
        self.BTN_Mantle_cambo4.addItem("Select a Model")
        self.BTN_Mantle_cambo4.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo4.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo4.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo4.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo4.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo4.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo4.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo4.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo4.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo4.addItem('Pyrolite5 (Ringwood, 1962*) ')  
        layout.addWidget(label1,i,0,1,1)  
        layout.addWidget(self.Depth4,i,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo4,i,2,1,2)

        i+=1
        label1 = QLabel("Layer5:")
        self.Depth5 = QLineEdit("")
        self.BTN_Mantle_cambo5=QComboBox()
        self.BTN_Mantle_cambo5.addItem("Select a Model")
        self.BTN_Mantle_cambo5.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo5.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo5.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo5.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo5.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo5.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo5.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo5.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo5.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo5.addItem('Pyrolite5 (Ringwood, 1962*) ')  
        layout.addWidget(label1,i,0,1,1)  
        layout.addWidget(self.Depth5,i,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo5,i,2,1,2)

        i+=1
        label1 = QLabel("Layer6:")
        self.Depth6 = QLineEdit("")
        self.BTN_Mantle_cambo6=QComboBox()
        self.BTN_Mantle_cambo6.addItem("Select a Model")
        self.BTN_Mantle_cambo6.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo6.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo6.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo6.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo6.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo6.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo6.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo6.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo6.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo6.addItem('Pyrolite5 (Ringwood, 1962*) ')  
        layout.addWidget(label1,i,0,1,1)  
        layout.addWidget(self.Depth6,i,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo6,i,2,1,2)

        
        self.setLayout(layout)
        self.setWindowTitle("Layers")
        self.resize(100,50)      
        
        
        self.okButton =QPushButton("&OK",self)
        self.okButton.clicked.connect(self.OK)
        self.okButton.setAutoDefault(False)
        self.cancelButton = QPushButton("Cancel",self)
        self.cancelButton.clicked.connect(self.Cancel)
        self.cancelButton.setAutoDefault(False)

        i+=1
        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.okButton)
        buttonLayout.addWidget(self.cancelButton)        
        layout.addLayout(buttonLayout, i, 0, 1, 2)
        
        self.Lines = [self.Depth1,self.Depth2,self.Depth3,self.Depth4,self.Depth5,self.Depth6]
        self.Combos = [self.BTN_Mantle_cambo1,self.BTN_Mantle_cambo2,self.BTN_Mantle_cambo3,self.BTN_Mantle_cambo4,self.BTN_Mantle_cambo5,self.BTN_Mantle_cambo6]        
        
    def OK(self):
        self.Change=True
        self.Depth=[]
        self.Phase=[]
        self.Phase_diagram=[]
        a=[]
        try:
            self.Depth.append(float(self.Depth1.text()));a.append(self.BTN_Mantle_cambo1.currentIndex())
        except:
            pass
        try:
            self.Depth.append(float(self.Depth2.text()));a.append(self.BTN_Mantle_cambo2.currentIndex())
        except:
            pass
        try:
            self.Depth.append(float(self.Depth3.text()));a.append(self.BTN_Mantle_cambo3.currentIndex())
        except:
            pass
        try:
            self.Depth.append(float(self.Depth4.text()));a.append(self.BTN_Mantle_cambo4.currentIndex())
        except:
            pass
        try:
            self.Depth.append(float(self.Depth5.text()));a.append(self.BTN_Mantle_cambo5.currentIndex())
        except:
            pass
        try:
            self.Depth.append(float(self.Depth6.text()));a.append(self.BTN_Mantle_cambo6.currentIndex())
        except:
            pass
        self.Result(a)
        self.CombosNumber=a
        for i in range(len(self.Depth)-1):
            if self.Depth[i+1] <= self.Depth[i] or self.Depth[i] >= 800:
                lalala
        if self.Depth[-1] <800 or self.Depth[0]<0:
            lalala
        if len(self.Depth) != len(self.Phase):
            lalala
        self.close()

    
    def Cancel(self):
        self.Change=False
        self.close()
        pass
        
    def Result(self,list1):
        for a in list1:
            if a == 1:
                Type = 'Basalt'
            if a == 2:
                Type = 'Harzburgite100'
            if a == 3:
                Type = 'Pyrolite'   
            if a == 4:
                Type = 'Piclogite'
            if a ==5:
                Type = 'Harzburgite90'       
            if a ==6:
                Type = 'Harzburgite80'          
            if a ==7:
                Type = 'Harzburgite70'
            if a ==8:
                Type = 'Harzburgite60'
            if a ==9:
                Type = 'Harzburgite50'
            if a ==10:
                Type = 'Harzburgite40' 
            if a ==11:
                Type = 'Harzburgite30'
            if a ==12:
                Type = 'Harzburgite20' 
            if a ==13:
                Type = 'Harzburgite10'
            if a ==14:
                Type = 'Pyrolite1'
            if a ==15:
                Type = 'Pyrolite2'
            if a ==16:
                Type = 'Pyrolite3'
            if a ==17:
                Type = 'Pyrolite4'
            if a ==18:
                Type = 'Pyrolite5'
            self.Phase.append(Type)


    def Reopen(self):
        if len(self.Depth) ==0:
            pass
        else:
            for i in range(len(self.Depth)):
                self.Lines[i].setText(str(self.Depth[i]))
                self.Combos[i].setCurrentIndex(self.CombosNumber[i])
                
    def Return(self):
        return self.Depth,self.Phase,self.CombosNumber                
        

if __name__ == "__main__":
    import sys
    layer_name=['Harzburgite10','Harzburgite100','Harzburgite10']
    layer_number=[1,1,1]
    layer_depth=[300,600,800]
    a = Layers(layer_name,layer_depth)
    
#==============================================================================
#     qApp = QApplication(sys.argv)
#     aw = GUI_layers(layer_depth,layer_name,layer_number)
#     aw.show()
#     sys.exit(qApp.exec_())       
#==============================================================================
