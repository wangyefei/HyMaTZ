# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 14:42:23 2016

@author: Fei
"""

import sys
import time
#from PyQt4 import QtGui,QtCore
try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import  *
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar    
except:
    from PyQt5.QtWidgets import (QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QHBoxLayout,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (QThread,pyqtSignal,PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar     
   

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.interpolate import interp1d



font = {'family' : 'serif',  
    'color'  : 'blue',  
    'weight' : 'normal',  
    'size'   : 10,  
    } 

   
    
def print_in_canvas_fit(fit_result,Harzburgite0,Harzburgite1):
    p,t,c,w1,w2,w3 = fit_result
    w1*=0.7;w2*=3.3; w3*=2.7
    print (p,t,c,w1,w2,w3)
    pyrolite_composition=dict()
    Mixture_composition = dict()
    Name=['CAO','AL2O3','NA2O','MGO','FEO','SIO2']
    for i in Name:
        pyrolite_composition[i]=Harzburgite0[i]*(1-0.82)+Harzburgite1[i]*0.82
        Mixture_composition[i]=Harzburgite0[i]*(1-0.82-c)+Harzburgite1[i]*(0.82+c)
    
    print (Name,(1-0.82-c))   
    string = 'Best fit:\n'+' delta Temperature (K):'+"{:5.2f}".format(t)+' water in Olivine (wt%)'+"{:5.2f}".format(w1)+' water in Wadsleyite (wt%)'+"{:5.2f}".format(w2)+' water in Ringwoodite (wt%)'+"{:5.2f}".format(w3)
    string1='Pyrolite composition: '+' CAO:'+"{:5.2f}".format(pyrolite_composition['CAO'])+' AL2O3:'+"{:5.2f}".format(pyrolite_composition['AL2O3'])+' NA2O:'+"{:5.2f}".format(pyrolite_composition['NA2O'])+' MGO:'+"{:5.2f}".format(pyrolite_composition['MGO'])+' FEO:'+"{:5.2f}".format(pyrolite_composition['FEO'])+' SIO2:'+"{:5.2f}".format(pyrolite_composition['SIO2'])  
    string2='new composition: '+' CAO:'+"{:5.2f}".format(Mixture_composition['CAO'])+' AL2O3:'+"{:5.2f}".format(Mixture_composition['AL2O3'])+' NA2O:'+"{:5.2f}".format(Mixture_composition['NA2O'])+' MGO:'+"{:5.2f}".format(Mixture_composition['MGO'])+' FEO:'+"{:5.2f}".format(Mixture_composition['FEO'])+' SIO2:'+"{:5.2f}".format(Mixture_composition['SIO2'])
    return string,string1,string2
    
    
def composition__output(composition,Harzburgite0=[ 13.88, 10.19,2.18,14.94,7.06,51.75],Harzburgite1=[0.81,0.53,0.,56.51,6.07,36.07]):
    pyrolite_composition=dict()
    Mixture_composition = dict()
    Name=['CAO','AL2O3','NA2O','MGO','FEO','SIO2']
    for i in range(6):
        pyrolite_composition[i]=Harzburgite0[i]*(1-0.82)+Harzburgite1[i]*0.82
        Mixture_composition[i]=Harzburgite0[i]*(1-composition/100)+Harzburgite1[i]*(composition/100)
      
    return pyrolite_composition,Mixture_composition
    
    
class PLOT(object):
    '''
    this plot
    '''
    def __init__(self,x,y,xlable='Depth (km)',ylable='y'):
        self.fig=plt.figure(dpi=200)
        ax2=self.fig.add_subplot(111)
        ax2.plot(x,y,color='r')
        #ax2.set_xlim(40,700)
        ax2.set_xlabel(xlable,fontdict=font)
        ax2.set_ylabel(ylable,fontdict=font)
        #ax2.legend(bbox_to_anchor=(0.05, 0.7), loc=3, borderaxespad=0.)
        
    def Back(self):
        return self.fig        


class wfQ_Canvas(FigureCanvas):
    """Class to represent the FigureCanvas widget"""
    def __init__(self, parent):
        # plot definition
        self.fig=Figure()
        FigureCanvas.__init__(self,self.fig)
        # set the parent widget
        self.setParent(parent)
        #self.plot()
        self.fig.canvas.draw()
        # we define the widget as expandable
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        # notify the system of updated policy
        FigureCanvas.updateGeometry(self)
              
  
        
class wfQ_ApplicationWindow(QDialog):
    """ main window"""
    def __init__(self,xx=None,x=[1,2,4],y=[2,3,4],name='None',xlabel='Depth (km)',ylabel='123'):
        
        ''' initialization of Qt MainWindow widget '''
        super(wfQ_ApplicationWindow, self).__init__()
        self.name=name
        self.xlabel=xlabel;self.ylabel=ylabel
        self.setWindowTitle(name)
        self.figure = Figure(dpi=200)
        self.canvas = FigureCanvas(self.figure)        
        self.ax = self.figure.add_subplot(111)        
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        self.buttonplot = QPushButton('plot')
        self.buttonplot.clicked.connect(self.plot)
        self.buttonimport = QPushButton('import')
        self.buttonimport.clicked.connect(self.import1)
        self.buttonexport = QPushButton('export')
        self.buttonexport.clicked.connect(self.export)
        self.buttonAdd = QPushButton('Add')
        self.buttonAdd.clicked.connect(self.Add)
        self.buttonMinus = QPushButton('Minus')
        self.buttonMinus.clicked.connect(self.Minus)

        self.AddEdit = QLineEdit('')
        self.AddEdit.returnPressed.connect(self.Add_num)
        #self.connect(self.AddEdit,QtCore.SIGNAL("returnPressed()"),)          
        
        self.MinusEdit = QLineEdit('')
        self.MinusEdit.returnPressed.connect(self.Minus_num)
        #self.connect(self.MinusEdit,QtCore.SIGNAL("returnPressed()"),)    
        
        layout = QGridLayout()
        layout.addWidget(self.canvas,1,1,5,6)
        layout.addWidget(self.AddEdit,7,1,1,3)        
        layout.addWidget(self.MinusEdit,7,4,1,3)
        layout.addWidget(self.buttonAdd,8,1,1,3)   
        layout.addWidget(self.buttonMinus,8,4,1,3)   
        layout.addWidget(self.toolbar,6,1,1,6)   
        layout.addWidget(self.buttonplot,9,1,1,6)
        layout.addWidget(self.buttonimport,10,1,1,6)
        layout.addWidget(self.buttonexport,11,1,1,6)
        
        self.setLayout(layout)
        
    def Add_num(self):
        self.add =  (float(self.AddEdit.text())) 
        for i in range(len(self.y)):
            self.y[i]+=self.add  
        self.plot()      
        
    def Minus_num(self):
        self.minus =  (float(self.MinusEdit.text())) 
        for i in range(len(self.y)):
            self.y[i]-=self.minus    
        self.plot()
            
    def change(self,x=[1,2,700],y=[2,3,4]):
        self.depth=x;self.y=y;

    def Add(self):
        pass

    def Minus(self):
        pass
        
    def export(self):
        fileName = QFileDialog.getSaveFileName(self)
        if fileName:
            file=open(fileName,'w')
            file.write('#  Depth');file.write('  ');file.write(self.name);file.write('\n')
            for i in range(len(self.depth)):
                file.write(str(self.depth[i]))
                file.write('  ');
                file.write(str(self.y[i]))
                file.write('\n')
        
    def import1 (self):
        fileName = QFileDialog.getOpenFileName(self)
        m=len(self.depth)
        file=open(fileName)
        print ('wf1')
        depth=[];y=[]
        for line in file:
            line=line.split()
            if line[0] == '#':
                pass
            else:
                depth.append(float(line[0]))
                y.append(float(line[1]))
        file.close()
        print ('wf2')
        depth=np.array(depth)
        y = np.array(y)
        print ('wf3')
        f_y = interp1d(depth,y)
        print ('wf4')
        self.y=f_y(self.depth)
        print ('wf5')
        self.plot()
            
    def output(self):
        return self.depth,self.y
        
    def plot(self):
        ''' plot some stuff '''
        self.ax = self.figure.add_subplot(111)
        self.ax.clear()
        self.ax.plot(self.depth,self.y)
        self.ax.set_xlabel(self.xlabel,fontdict=font)
        self.ax.set_ylabel(self.name,fontdict=font)
        self.canvas.draw()   

        

class wfQ_Phase_plot(QWidget):
    """ main window"""
    def __init__(self,xx=None,fig=None):
        
        ''' initialization of Qt MainWindow widget '''
        super(wfQ_Phase_plot, self).__init__()
        self.figure = fig
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
 
        self.button = QPushButton('Plot')
        #self.button.clicked.connect(self.plot)

        layout = QGridLayout()
        layout.addWidget(self.canvas,1,1,5,5)
        layout.addWidget(self.toolbar,6,1,1,5)       
        layout.addWidget(self.button,7,1,1,5)
        
        self.setLayout(layout)
        
class Test_Thread(QThread):
    finished = pyqtSignal()

    def __init__(self):
        QThread.__init__(self)

    def run(self):
        time.sleep(10)

        self.finished.emit()
        self.terminate()  
    
    def Finished(self):
        self.finished.emit()
        self.terminate()         



class WFQ_condition(QWidget):
    
    def __init__(self,D,P,T,C,W):
        super(WFQ_condition, self).__init__()
        self.figure = plt.figure(dpi=200)
        self.canvas = FigureCanvas(self.figure)
        self.ax1 = self.figure.add_subplot(411) 
        self.ax1.plot(D,P)
        self.ax2 = self.figure.add_subplot(412) 
        self.ax2.plot(D,T)
        self.ax3 = self.figure.add_subplot(413) 
        self.ax3.plot(D,C)
        self.ax4 = self.figure.add_subplot(414) 
        self.ax4.plot(D,W)        
        self.toolbar = NavigationToolbar(self.canvas, self)
 
        self.button = QPushButton('Plot')


        layout = QGridLayout()
        layout.addWidget(self.canvas,1,1,5,5)
        layout.addWidget(self.toolbar,6,1,1,5)       
        layout.addWidget(self.button,7,1,1,5)
        self.setLayout(layout)
             
class MantleDlg(QDialog):

    def __init__(self, choice1=None,choice2=None,precentage=None):
        super(MantleDlg, self).__init__()

        Label1 = QLabel("Endmember1:")
        Label2 = QLabel("Endmember2:")
        
        #self.choice1 = choice1
        #self.choice2 = choice2
        #self.precentage = precentage
        
        
        self.Harzburgite = np.array([0.81, 0.53, 0.00, 56.51, 6.07, 36.07])
        self.Basalt      = np.array([13.88, 10.19, 2.18, 14.94, 7.06, 51.75])
        self.Pyrolite    = np.array([2.93553,2.21491,0.109863,49.8509,6.16912,38.6950])
        
        self.BTN_Mantle_cambo=QComboBox()
        self.BTN_Mantle_cambo.addItem("Select a Model")
        self.BTN_Mantle_cambo.addItem("Basalt(Xu et al. 2008)")
        self.BTN_Mantle_cambo.addItem("Harzburgite (Xu et al. 2008)")
        self.BTN_Mantle_cambo.addItem("Pyrolite (Xu et al. 2008)")
        self.BTN_Mantle_cambo.addItem('Piclogite (Weidner, 1986) ')
        self.BTN_Mantle_cambo.addItem("90%molHarzburgite + 10%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("80%molHarzburgite + 20%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("70%molHarzburgite + 30%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("60%molHarzburgite + 40%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("50%molHarzburgite + 50%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("40%molHarzburgite + 60%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("30%molHarzburgite + 70%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("20%molHarzburgite + 80%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("10%molHarzburgite + 90%molBasalt (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem('Pyrolite1 (McDonough and Sun, 1995) ')
        self.BTN_Mantle_cambo.addItem('Pyrolite2 (Jagoutz et al., 1979) ')
        self.BTN_Mantle_cambo.addItem('Pyrolite3 (Green et al., 1979) ')
        self.BTN_Mantle_cambo.addItem('Pyrolite4 (Taylor and McLennan, 1985) ')
        self.BTN_Mantle_cambo.addItem('Pyrolite5 (Ringwood, 1962*) ')       
        if choice1 is not None:
            self.BTN_Mantle_cambo.setCurrentIndex(choice1)
        else:
            self.BTN_Mantle_cambo.setCurrentIndex(self.BTN_Mantle_cambo.findText("Select a Model"))
        self.BTN_Mantle_cambo.currentIndexChanged.connect(self.Model_select)
        
        
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
        if choice1 is not None:
           self.BTN_Mantle_cambo1.setCurrentIndex(choice2)
        else:
            self.BTN_Mantle_cambo1.setCurrentIndex(0)
        self.BTN_Mantle_cambo1.currentIndexChanged.connect(self.Model_select1)


        self.Label = QLabel("&Endmember2 precentage (volume):")
        self.Mixture = QSpinBox()
        self.Label.setBuddy(self.Mixture)
        self.Mixture.valueChanged.connect(self.valuechange)
        self.Mixture.setMaximum(100)
        self.Mixture.setMinimum(0)
        try:
            self.Mixture.setValue(precentage)
        except:
            pass
        
        self.okButton =QPushButton("&OK",self)
        self.okButton.clicked.connect(self.OK)
        self.cancelButton = QPushButton("Cancel",self)
        self.cancelButton.clicked.connect(self.Cancel)

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.okButton)
        buttonLayout.addWidget(self.cancelButton)
        
        layout = QGridLayout()
        layout.addWidget(Label1,0,0,1,1)
        layout.addWidget(Label2,1,0,1,1)
        layout.addWidget(self.Label,2,0,1,1)
        layout.addWidget(self.Mixture,2,1,1,1)
        layout.addLayout(buttonLayout, 3, 0, 1, 2)
        layout.addWidget(self.BTN_Mantle_cambo,1,1,1,1)
        layout.addWidget(self.BTN_Mantle_cambo1,0,1,1,1)
        self.setLayout(layout)
        self.setWindowTitle("Mantle Properties")
        self.resize(100,50)

    
    def valuechange(self):
        self.value=self.Mixture.value()
        #self.Label.setText("Hargurbite precentage (%):") 
    def Model_select(self,a):
        if a == 1:
            self.Type = 'Basalt'
            self.Model_composition = self.Basalt
        if a == 2:
            self.Type = 'Harzburgite100'
            self.Model_composition = self.Harzburgite
        if a == 3:
            self.Type = 'Pyrolite'   
            self.Model_composition =  self.Pyrolite
        if a == 4:
            self.Type = 'Piclogite'
            self.Model_composition = np.array([  7.22,8.16,2.80,30.14,4.58,47.07])
        if a ==5:
            self.Type = 'Harzburgite90' 
            Harzburgite_precentage = 0.9
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt               
        if a ==6:
            self.Type = 'Harzburgite80'
            Harzburgite_precentage = 0.8
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt              
        if a ==7:
            self.Type = 'Harzburgite70'
            Harzburgite_precentage = 0.7
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==8:
            self.Type = 'Harzburgite60'
            Harzburgite_precentage = 0.6
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==9:
            self.Type = 'Harzburgite50'
            Harzburgite_precentage = 0.5
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==10:
            self.Type = 'Harzburgite40' 
            Harzburgite_precentage = 0.4
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==11:
            self.Type = 'Harzburgite30'
            Harzburgite_precentage = 0.3
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==12:
            self.Type = 'Harzburgite20' 
            Harzburgite_precentage = 0.2
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==13:
            self.Type = 'Harzburgite10'
            Harzburgite_precentage = 0.1
            self.Model_composition = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==14:
            self.Type = 'Pyrolite1'
            self.Model_composition = np.array([  3.31,2.28,0.30,49.06,5.86,39.18])
        if a ==15:
            self.Type = 'Pyrolite2'
            self.Model_composition = np.array([  3.25,2.03,0.27,49.56,5.67,39.19])
        if a ==16:
            self.Type = 'Pyrolite3'
            self.Model_composition = np.array([  3.14,2.23,0.33,49.94,5.48,38.85])
        if a ==17:
            self.Type = 'Pyrolite4'
            self.Model_composition = np.array([  2.71,1.87,0.28, 45.73,5.83,43.55])
        if a ==18:
            self.Type = 'Pyrolite5'
            self.Model_composition = np.array([  2.90,1.82,0.50, 48.96,6.22,39.58])
        self.a1 = a
        
    def Model_select1(self,a):
        if a == 1:
            self.Type1 = 'Basalt'
            self.Model_composition1 = self.Basalt
        if a == 2:
            self.Type1 = 'Harzburgite100'
            self.Model_composition1 = self.Harzburgite
        if a == 3:
            self.Type1 = 'Pyrolite'   
            self.Model_composition1 =  self.Pyrolite
        if a == 4:
            self.Type1 = 'Piclogite'
            self.Model_composition1 = np.array([  7.22,8.16,2.80,30.14,4.58,47.07])
        if a ==5:
            self.Type1 = 'Harzburgite90' 
            Harzburgite_precentage = 0.9
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt               
        if a ==6:
            self.Type1 = 'Harzburgite80'
            Harzburgite_precentage = 0.8
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt              
        if a ==7:
            self.Type1 = 'Harzburgite70'
            Harzburgite_precentage = 0.7
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==8:
            self.Type1 = 'Harzburgite60'
            Harzburgite_precentage = 0.6
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==9:
            self.Type1 = 'Harzburgite50'
            Harzburgite_precentage = 0.5
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==10:
            self.Type1 = 'Harzburgite40' 
            Harzburgite_precentage = 0.4
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==11:
            self.Type1 = 'Harzburgite30'
            Harzburgite_precentage = 0.3
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==12:
            self.Type1 = 'Harzburgite20' 
            Harzburgite_precentage = 0.2
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==13:
            self.Type1 = 'Harzburgite10'
            Harzburgite_precentage = 0.1
            self.Model_composition1 = Harzburgite_precentage*self.Harzburgite + (1-Harzburgite_precentage)*self.Basalt
        if a ==14:
            self.Type1 = 'Pyrolite1'
            self.Model_composition1 = np.array([  3.31,2.28,0.30,49.06,5.86,39.18])
        if a ==15:
            self.Type1 = 'Pyrolite2'
            self.Model_composition1 = np.array([  3.25,2.03,0.27,49.56,5.67,39.19])
        if a ==16:
            self.Type1 = 'Pyrolite3'
            self.Model_composition1 = np.array([  3.14,2.23,0.33,49.94,5.48,38.85])
        if a ==17:
            self.Type1 = 'Pyrolite4'
            self.Model_composition = np.array([  2.71,1.87,0.28, 45.73,5.83,43.55])
        if a ==18:
            self.Type = 'Pyrolite5'
            self.Model_composition1 = np.array([  2.90,1.82,0.50, 48.96,6.22,39.58])
        self.a2=a
        
    def OK(self):
        self.composition = self.Model_composition*self.value/100+self.Model_composition1*(1-self.value/100)
        self.close()
        
    def Cancel(self):
        self.close()
        
        

class WfQ_profile_control(QWidget):

    def __init__(self,Depth=[100,400,600],control1=[0,1,2],control2=[0,2,3],control3=[1,4,5]\
                                         ,control11=[1,0,0],control22=[0,1,0],control33=[0,0,1]):
        super(WfQ_profile_control,self).__init__()
        
        self.Depth = Depth
        self.control1 = np.array(control1)
        self.control2 = np.array(control2)
        self.control3 = np.array(control3)
        self.control11 = np.array(control11)
        self.control22 = np.array(control22)
        self.control33 = np.array(control33)
        
        self.figure = plt.figure(dpi=200)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        layout = QGridLayout()
        layout.addWidget(self.canvas,1,1,5,5)
        layout.addWidget(self.toolbar,6,1,1,5)  

        Label= QLabel()
        Label.setText('OL water capacity') 
        layout.addWidget(Label,7,1,1,2)
        Label= QLabel()
        Label.setText('WA water capacity') 
        layout.addWidget(Label,7,3,1,2)
        Label= QLabel()
        Label.setText('RI water capacity') 
        layout.addWidget(Label,7,5,1,2)
        
        self.textEdit_OL = QLineEdit()
        self.textEdit_OL.setFixedWidth(200)
        self.textEdit_OL.setText("")
        self.textEdit_OL.returnPressed.connect(self.editor_OL)
        layout.addWidget(self.textEdit_OL,8,1,1,2)
        #self.connect(self.textEdit_OL,QtCore.SIGNAL("returnPressed()"),)
        self.textEdit_WA = QLineEdit()
        self.textEdit_WA.setFixedWidth(200)
        self.textEdit_WA.setText("")
        self.textEdit_WA.returnPressed.connect(self.editor_WA)
        layout.addWidget(self.textEdit_WA,8,3,1,2)
        #self.connect(self.textEdit_WA,QtCore.SIGNAL("returnPressed()"),self.editor_WA)        
        self.textEdit_RI = QLineEdit()
        self.textEdit_RI.setFixedWidth(200)
        self.textEdit_RI.setText("")
        self.textEdit_RI.returnPressed.connect(self.editor_RI)
        layout.addWidget(self.textEdit_RI,8,5,1,2)
        #self.connect(self.textEdit_RI,QtCore.SIGNAL("returnPressed()"),self.editor_RI)
        
 
        self.button = QPushButton('Plot')
        self.button.clicked.connect(self.Plot)    
        layout.addWidget(self.button,9,1,1,5)       
        self.setLayout(layout)
        
    def editor_OL(self):
        self.OL=(float(self.textEdit_OL.text())) 
        for i in range(len(self.control1)):
            if self.control11[i] != 0:
                self.control1[i] = self.OL
        pass
    def editor_WA(self):
        self.WA=(float(self.textEdit_WA.text())) 
        for i in range(len(self.control1)):
            if self.control22[i] != 0:
                self.control2[i] = self.WA  
        pass
    def editor_RI(self):
        self.RI=(float(self.textEdit_RI.text())) 
        for i in range(len(self.control1)):
            if self.control33[i] != 0:
                self.control3[i] = self.RI     
        pass   
    
    def Plot(self):
        self.ax.clear()
        self.ax.plot(self.Depth,self.control1,'r')
        self.ax.plot(self.Depth,self.control2,'b')
        self.ax.plot(self.Depth,self.control3,'g')
        #self.ax.set_xlabel(self.xlabel,fontdict=font)
        self.ax.set_ylabel('water wt%',fontdict=font)
        #self.ax.le
        self.ax.set_ylim(0,4)
        self.canvas.draw()  
    



        
        

    
if __name__ == "__main__":   
    qApp = QApplication(sys.argv)
    # instantiate the ApplicationWindow widget
    aw = MantleDlg(choice1=4,choice2=2,precentage=11)
    #aw.change()
    # show the widget
    aw.show()
    
    # start the Qt main loop execution, exiting from this script
    # with the same return code of Qt application
    sys.exit(qApp.exec_())       
   