# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 15:49:49 2017

@author: Fei
"""

import os
import platform
import sys
try:
    from PyQt4.QtCore import (PYQT_VERSION_STR, QFile, QFileInfo, QSettings,
            QT_VERSION_STR, QTimer, QVariant, Qt, SIGNAL)
    PYQT= 4
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
    PYQT= 5
    from PyQt5.QtCore import (PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar     
    
#from GUI_tools.wf_GUI_plot import *
#from Mineral_Physics.Rock import *
from Mineral_Physics.Stix2011data import *

#SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
#SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

class Mineral_Params(QDialog,QTableView):
    
    def __init__(self,Pressure=1e5,Temperature=300,Volume=None,VolumeOff=True):
        super(Mineral_Params,self).__init__()
        self.Solidsolution = [an,ab,sp,hc,en,fs,mgts,odi,hpcen,hpcfs,di,he,cen,cats,jd,py,al,gr,mgmj,jdmj, \
                          capv,fo,fa,mgwa,fewa,mgri,feri,mgil,feil,co,mgpv,fepv,alpv,mppv,fppv,appv,mgcf,fecf, \
                          nacf,pe,wu,qtz,coes,st,an,ky,neph]
        self.row=len(self.Solidsolution)
        self.layout = QGridLayout(self) 
        
        self.Pressure=Pressure;self.Temperature=Temperature

        self.table = QTableView()
        self.model = QStandardItemModel(10,10,self)
        if Volume is None:
            self.Volume = Volume
            self.model.setHorizontalHeaderLabels(['Name','formuala','Vp (km/s)','Vs (km/s)','Rho (kg/m³)','F (KJ/mol)','V (cm³/mol)','K (GPa)',"k'",'θ'+'(K)','γ','q','G (GPa)',"G'",'η'])
        else:
            self.Volume=Volume[:-3]
            self.model.setHorizontalHeaderLabels(['Name','formuala','Vp (km/s)','Vs (km/s)','Rho (kg/m³)','F (KJ/mol)','V (cm³/mol)','K (GPa)',"k'",'θ'+'(K)','γ','q','G (GPa)',"G'",'η','Volume fraction'])
        

        #self.model.setHorizontalHeaderLabels(['Name','formuala','Vp (km/s)','Vs (km/s)','Rho (kg/m³)','F (KJ/mol)','V (cm³/mol)','K (GPa)',"k'",'θ'+'(K)','γ','q','G (GPa)',"G'",'η','Volume fraction'])
            
        self.table.setModel(self.model)
                
        self.layout.addWidget(self.table,0,0,1,5) 
        
        self.BTN()
        self.table_view()
        self.resize(1800, 800)
        #self.restore()

        
    def BTN(self):
        
        self.Update = QPushButton(self)
        self.Update.setText("Update")
        self.Update.clicked.connect(self.update_data)
        self.Update.setAutoDefault(False)
                     
        self.Restore = QPushButton(self)
        self.Restore.setText("Restore")
        self.Restore.clicked.connect(self.restore)  
        self.Restore.setAutoDefault(False)
        
        self.Random = QPushButton(self)
        self.Random.setText("Random within error")
        self.Random.clicked.connect(self.random)
        self.Random.setAutoDefault(False)

        self.Save = QPushButton(self)
        self.Save.setText("Save")
        self.Save.clicked.connect(self.save)
        self.Save.setAutoDefault(False)
        
        self.Iutput = QPushButton(self)
        self.Iutput.setText("Load")
        self.Iutput.clicked.connect(self.iutput)
        self.Iutput.setAutoDefault(False)
        
        self.layout.addWidget(self.Update,1,0,1,1)
        self.layout.addWidget(self.Restore,1,1,1,1)
        self.layout.addWidget(self.Random,1,2,1,1)
        self.layout.addWidget(self.Save,1,3,1,1)
        self.layout.addWidget(self.Iutput,1,4,1,1)

    def save(self):
        options = QFileDialog.Options()
        fileName = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)   
        #except:
        #    fileName, _= QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)   
        if PYQT == 5:
            fileName = fileName[0]
        else:
            pass
        if fileName:
            file1 = open(fileName,'w')
            #print ('open',fileName)
            file1.write('name, F_0,V_0,K_0,Kprime_0 Debye_0,grueneisen-0,q_0,G_0,Gprime_0, eta_s0 ' )
            file1.write('\n')
            for i in self.Solidsolution:
                for j in i.parameters(self.Pressure,self.Temperature):
                    file1.write('{:8.8}'.format(j))
                    file1.write(' ')
                file1.write('\n')
            file1.close()
            #print ('save')
        
    def iutput(self):
        fileName = QFileDialog.getOpenFileName(self, 'Open file', '/home') 
        if PYQT == 5:
            fileName = fileName[0]
        else:
            pass
        if fileName:
           file1 = open(fileName,'r') 
           next(file1)
           for i,line in enumerate(file1):
                line = line.split()
                newItem = QStandardItem(self.Solidsolution[i].name) 
                self.model.setItem(i, 0, newItem) 
                for col in range(len(line)):
                    try:
                        string = "{:5.2f}".format(line[col])
                    except:
                        string = line[col]
                    item = QStandardItem(string)
                    self.model.setItem(i, col, item) 
           file1.close()
                    #print (string)
        #print ('input')
                
        #self.table_view()               
            
    def table_view(self):
        #self.row_list=[]
        for i in range(self.row):
            newItem = QStandardItem(self.Solidsolution[i].name) 
            self.model.setItem(i, 0, newItem) 
            params=self.Solidsolution[i].parameters(1e5,300)
            for col in range(len(params)):
                try:
                    string = "{:5.2f}".format(params[col])
                except:
                    string = params[col]
                item = QStandardItem(string)
                self.model.setItem(i, col, item)
        if self.Volume is not None:
            for num,i in enumerate(self.Volume):
                newItem = QStandardItem(str(self.Volume[num]))
                self.model.setItem(num, 15, newItem) 
                
    def random(self):
        for i in range(self.row):
            newItem = QStandardItem(self.Solidsolution[i].name) 
            self.model.setItem(i, 0, newItem) 
            params=self.Solidsolution[i].parameters_random(self.Pressure,self.Temperature)
            #print ('random')
            for col in range(len(params)):
                try:
                    string = "{:5.2f}".format(params[col])
                except:
                    string = params[col]
                item = QStandardItem(string)
                self.model.setItem(i, col, item)     
            
        #self.table_view()
        #print ('random')
        
    def update_data(self):
        for i in range(len(self.Solidsolution)):
            params=[]
            for j in range(15):
                index=self.model.index(i,j)
                params.append(self.model.itemData(index)[0])
            self.Solidsolution[i].change_parameters(params)
            a,b,c=self.Solidsolution[i].Vp_Vs(self.Pressure*1e4,self.Temperature);c*=1000
            item=QStandardItem("{:5.2f}".format(a));self.model.setItem(i, 2, item)
            item=QStandardItem("{:5.2f}".format(b));self.model.setItem(i, 3, item)
            item=QStandardItem("{:5.2f}".format(c));self.model.setItem(i, 4, item)        
        #print ('wf')
        
    def restore(self):
        try:
            address=os.path.join(os.path.dirname(__file__),'EXPDATA','All.txt')
            self.address=address
            file1=open(address,'r+')
        except:
            address=os.path.join(os.path.dirname(__file__),'Mineral_Physics','EXPDATA','All.txt')
            self.address=address
            file1=open(address,'r') 
            next(file1)
            for i,line in enumerate(file1):
                line = line.split()
                newItem = QStandardItem(self.Solidsolution[i].name) 
                self.model.setItem(i, 0, newItem) 
                for col in range(len(line)):
                    try:
                        string = "{:5.2f}".format(line[col])
                    except:
                        string = line[col]
                    item = QStandardItem(string)
                    self.model.setItem(i, col, item)   
            file1.close()
        self.update_data()
        #print ('wf')

        
                
                
if __name__ == "__main__":
     qApp = QApplication(sys.argv)
     aaa=[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   4.84480000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         7.72000000e-02,   2.02300000e-01,   4.79200000e-01,
         5.66080000e+01,   2.90520000e+00,   7.87200000e-01,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         1.05780000e+00,   4.08000000e-02,   2.63200000e-01,
         2.65540000e+01,   6.17790000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00]
     aw = Mineral_Params(1e5,300,VolumeOff=True)
     aw.show()
     sys.exit(qApp.exec_())