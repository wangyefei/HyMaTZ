# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 10:06:34 2017

@author: Fei
"""

import os
import platform
import sys
import copy
import numpy as np
try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar      
except:
    from PyQt5.QtWidgets import (QAbstractItemView,QVBoxLayout,QSizePolicy,QProgressDialog,QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QHBoxLayout,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (QThread,pyqtSignal,PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar  

import  matplotlib.pyplot  as plt 
from matplotlib.figure import Figure   
from GUI_tools.wf_GUI_plot import (Test_Thread,wfQ_ApplicationWindow,MantleDlg)
from Mineral_Physics.EOS_modify import EOS,test

class GUI_table(QDialog):
    
    def __init__(self,params=test.HPparams ):
        super(GUI_table,self).__init__()
        
        self.Main = QGridLayout(self)
        self.resize(1600,400)
        self.Main_Splitter = QSplitter(Qt.Vertical)

        self.table = QTableView()
        
        self.header = list(params.keys())
        self.value = list(params.values())
        self.model = QStandardItemModel(1,len(self.header),self)
        self.model.setHorizontalHeaderLabels(self.header)
        for i in range(len(self.header)):
                newItem = QStandardItem(str(self.value[i])) 
                self.model.setItem(0, i, newItem) 
                
        self.table.setModel(self.model)   

        
        self.Main_Splitter.addWidget(self.table)
        self.Main.addWidget(self.Main_Splitter,0,0,2,16)

        self.okButton =QPushButton("&OK",self)
        self.okButton.clicked.connect(self.OK)
        self.cancelButton = QPushButton("Cancel",self)
        self.cancelButton.clicked.connect(self.Cancel)

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.okButton)
        buttonLayout.addWidget(self.cancelButton)
        self.Main.addLayout(buttonLayout,2,14,1,2)        
 

    def OK(self):
        self.get_table_data()
        self.close()
        
    def Cancel(self):
        self.close()


    def get_table_data(self):              
        self.params=np.zeros(len(self.header))
        for j in range(len(self.header)):
            try:
                index=self.model.index(0,j)
                self.params[j]=(self.model.itemData(index)[0])        
            except:
                pass
                #a = self.showdialog('Please fill in all the blank ')
        return self.params


    def showdialog(self,message1=' ',message2=' '):
       msg = QMessageBox()
       msg.setIcon(QMessageBox.Information)    
       msg.setText(message1)
       msg.setWindowTitle("MessageBox")
       msg.setDetailedText(message2)
       msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
       if msg.exec_()== QMessageBox.Yes:
           return True
       else:
           return False


class GUI_user_input(QDialog):
    
    def __init__(self,D=np.array([1,2,3]),P=np.array([1,2,3]),T=np.array([1,2,3])):
        super(GUI_user_input, self).__init__() 
        self.D=np.array(D);self.P=np.array(P);self.T=np.array(T)
        self.M1=test
        self.inputnumber=0
        
        self.Main = QHBoxLayout(self)
        self.Main_Splitter = QSplitter(Qt.Horizontal)

        
        self.Plot_layout()
        self.Control_layout()
        
        self.Main_Splitter.addWidget(self.Control_Widget)
        self.Main_Splitter.addWidget(self.Plot_Widget)
        self.Main.addWidget(self.Main_Splitter)  

         
    def Control_layout(self):
        i=0
        self.Control_Widget = QWidget()
        self.BTN_layout = QGridLayout(self.Control_Widget)
        
        label1 = QLabel()
        label1.setText('1. Select a mineral or add a new one')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)
        
        i+=1
        self.Mineral_list=['Choose one','add',test]
        self.Mineral_list_Combox= QComboBox()
        self.Mineral_list_Combox.addItem('Choose one')
        self.Mineral_list_Combox.addItem('add')
        self.Mineral_list_Combox.addItem('example')
        self.Mineral_list_Combox.setCurrentIndex(self.Mineral_list_Combox.findText("Choose one"))
        self.Mineral_list_Combox.currentIndexChanged.connect(self.AddMineralList)
        self.BTN_layout.addWidget(self.Mineral_list_Combox,i,0,1,2)
        
        i+=1
        label1 = QLabel()
        label1.setText('2. Enter new name or change name')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)    
        
        i+=1
        self.EditName = QLineEdit()
        try:
            self.EditName.setText(self.MineralName)
        except:
            self.EditName.setText('')
        self.EditName.setStyleSheet("background-color: rgb(255, 255, 255)")    
        self.EditName.returnPressed.connect(self.editName)
        self.BTN_layout.addWidget(self.EditName,i,0,1,2)        
        
        i+=1
        label1 = QLabel()
        label1.setText('3. Enter new formula or change formula')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)          

        i+=1
        self.EditFormula = QLineEdit()
        try:
            self.EditFormula.setText(self.MineralFormula)
        except:
            self.EditFormula.setText('')
        self.EditFormula.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.EditFormula.returnPressed.connect(self.editFormula)
        self.BTN_layout.addWidget(self.EditFormula,i,0,1,2)   
        
        i+=1
        label1 = QLabel()
        label1.setText('4. Choose an EOS')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)    
        
        i+=1
        self.EOS_list=['Choose one','BM-MGD (Stixrude & Stixrude & Lithgow-Bertelloni (2005)','TEOS Holand & Powell (2011)','BM2rd','BM3rd']
        self.EOS_list_Combox= QComboBox()
        self.EOS_list_Combox.addItems(self.EOS_list)
        self.EOS_list_Combox.setCurrentIndex(self.EOS_list_Combox.findText("Choose one"))
        self.EOS_list_Combox.currentIndexChanged.connect(self.eos_choice) 
        self.BTN_layout.addWidget(self.EOS_list_Combox,i,0,1,2) 
        
        i+=1
        label1 = QLabel()
        label1.setText('5. Enter or change parameters')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)


        i+=1
        self.Parameter_Button1=QPushButton()
        self.Parameter_Button1.setText('Enter/Modify thermoelastic data')
        self.Parameter_Button1.clicked.connect(self.Viewthermoelasticdata)  
        self.Parameter_Button1.setAutoDefault(False)
        self.BTN_layout.addWidget(self.Parameter_Button1,i,0,1,2)
        
        i+=1
        self.Button1=QPushButton()
        self.Button1.setText('Comfirm input')
        self.Button1.clicked.connect(self.Comfirm)  
        self.Button1.setAutoDefault(False)
        self.BTN_layout.addWidget(self.Button1,i,0,1,2)        

        i+=1
        label1 = QLabel()
        label1.setText('7. Change pressure and temperature')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)

        i+=1
        self.Mantle_Properties_Combobox = QComboBox()
        self.Mantle_Properties_Combobox.addItem('Choose one')
        self.Mantle_Properties_Combobox.addItem('Pressure')
        self.Mantle_Properties_Combobox.addItem('Temperature')
        self.Mantle_Properties_Combobox.setCurrentIndex(self.Mantle_Properties_Combobox.findText("Choose one"))
        self.Mantle_Properties_Combobox.currentIndexChanged.connect(self.mantle_Properties_Button) 
        self.BTN_layout.addWidget(self.Mantle_Properties_Combobox,i,0,1,2) 
        
        i+=1
        label1 = QLabel()
        label1.setText('8. Make some plots')
        label1.setWordWrap(True)
        self.BTN_layout.addWidget(label1,i,0,1,2)  
        
        i+=1
        self.Plot_Combobox=QComboBox()
        self.PLOT=['Choose one','clean','Vp','Vs','Vp+Vs','Rho','K','G']
        self.Plot_Combobox.addItems(self.PLOT)
        self.Plot_Combobox.setCurrentIndex(self.Plot_Combobox.findText("Choose one"))
        self.Plot_Combobox.currentIndexChanged.connect(self.Plot_option)
        self.BTN_layout.addWidget(self.Plot_Combobox,i,0,1,2) 
        
        i+=1
        self.PLOT_Button=QPushButton()
        self.PLOT_Button.setText('Plot')
        self.PLOT_Button.clicked.connect(self.PLOT_function)
        self.PLOT_Button.setAutoDefault(False)        
        self.BTN_layout.addWidget(self.PLOT_Button,i,0,1,2) 

    def Change_MineralList(self,MineralList):
        self.Mineral_list = MineralList
        print ('1')
        self.Mineral_list_Combox= QComboBox()
        print ('2')
        self.Mineral_list_Combox.addItem('Choose one')
        self.Mineral_list_Combox.addItem('add')
        self.Mineral_list_Combox.addItem('example')
        print ('3')
        for i in self.Mineral_list[3:]:
            self.Mineral_list_Combox.addItem(i.name)
        print ('4')
        self.Mineral_list_Combox.setCurrentIndex(self.Mineral_list_Combox.findText("Choose one"))   
        print ('5')
        self.Mineral_list_Combox.currentIndexChanged.connect(self.AddMineralList)
        self.BTN_layout.addWidget(self.Mineral_list_Combox,1,0,1,2)
        print ('6')

    def AddMineralList(self,a):
        self.EditName.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.EditFormula.setStyleSheet("background-color: rgb(255, 255, 255)")
        
        if a == 0:
            self.inputnumber=a
            pass
        if a ==1:
            self.Mineral_list.append(EOS())   
            self.EditName.setText(str('input name'))
            self.EditFormula.setText(str('input formula'))
            self.inputnumber = a
            
        if a!=1 and a!=0:
            self.inputnumber=a
            self.MineralName = self.Mineral_list[a].name
            self.MineralFormula = self.Mineral_list[a].formula    
            self.EditName.setText(self.MineralName)
            self.EditFormula.setText(self.MineralFormula)

            
            
        #self.Mineral_list_Combox.setCurrentIndex(self.Mineral_list_Combox.findText("Choose one"))


    def editName(self):
        self.MineralName = ((self.EditName.text())) 
        if self.inputnumber >1:
            self.Mineral_list[self.inputnumber].Change_Name(self.MineralName)
            self.EditName.setStyleSheet("background-color: rgb(155, 155, 155)")
        if self.inputnumber ==1:
            self.Mineral_list[-1].Change_Name(self.MineralName)
            self.EditName.setStyleSheet("background-color: rgb(155, 155, 155)")            
            
    
    def editFormula(self):
        self.MineralFormula = ((self.EditFormula.text())) 
        if self.inputnumber >1:
            self.Mineral_list[self.inputnumber].Change_formula(self.MineralFormula)
            self.EditFormula.setStyleSheet("background-color: rgb(155, 155, 155)")
        if self.inputnumber ==1:
            self.Mineral_list[-1].Change_formula(self.MineralFormula)
            self.EditFormula.setStyleSheet("background-color: rgb(155, 155, 155)")          

         
    def eos_choice(self,a):
        self.eos=a
        pass
    
    def Viewthermoelasticdata(self):
        if self.inputnumber !=0 and self.inputnumber !=1:
            aa=self.Mineral_list[self.inputnumber]
        else:
            aa=test
        if self.eos == 1:
            self.viewthermo = GUI_table(aa.SLBparams)
        if self.eos == 2:
            self.viewthermo = GUI_table(aa.HPparams)    
        if self.eos == 3:
            self.viewthermo = GUI_table(aa.BM2rdparams)                    
        if self.eos == 4:
            self.viewthermo = GUI_table(aa.BM3rdparams)  
        if self.eos == 5:
            self.viewthermo = GUI_table(aa.BMparams)  
        
        print ('wf1')
        if self.viewthermo.exec_():
            pass
        
        self.result = copy.deepcopy(self.viewthermo.params)

        if self.eos == 1:
            self.Mineral_list[self.inputnumber].Change_SLB(self.result)
        if self.eos == 2:
            self.Mineral_list[self.inputnumber].Change_HP(self.result)  
        if self.eos == 3:
            self.Mineral_list[self.inputnumber].Change_BM2rd(self.result)                 
        if self.eos == 4:
            self.Mineral_list[self.inputnumber].Change_BM3rd(self.result)
        if self.eos == 5:
            self.Mineral_list[self.inputnumber].Change_BM(self.result)
            
    
    def Comfirm(self):
        if self.inputnumber >1:
            self.Mineral_list[self.inputnumber].Change_Name_formula(str(self.MineralName),str(self.MineralFormula))
        if self.inputnumber == 1:
            self.editName()
            self.editFormula()
            self.Mineral_list[-1].Change_Name_formula(str(self.MineralName),str(self.MineralFormula))
            self.Mineral_list_Combox.addItem(str(self.MineralName))
            self.Mineral_list_Combox.setCurrentIndex(self.Mineral_list_Combox.findText(str(self.MineralName)))
        print ('comfirm')    
    
    def mantle_Properties_Button(self,a):
        if a==1:
            self.Pressure= wfQ_ApplicationWindow(name="Pressure (Bar)")      
            self.Pressure.change(x=self.D,y=self.P)
            if self.Pressure.exec_():   
                self.D,self.P=self.Pressure.output()
        if a==2:
            self.Temperature= wfQ_ApplicationWindow(name="Temperature (K)")                     
            self.Temperature.change(x=self.D,y=self.T)
            if self.Temperature.exec_():  
                self.D,self.T=self.Temperature.output()       
        self.Mantle_Properties_Button.setCurrentIndex(self.Mantle_Properties_Button.findText("Choose one"))


   
    
    def Plot_option(self,a):
        self.plotoption=a
        
    
    def PLOT_function(self):
        Vp=[];Vs=[];Rho=[];K=[];G=[]
        if self.plotoption ==1:
            self.ax.cla()
            
        if self.eos==1:
            for i in range(len(self.D)):
                vp,vs,rho,k,g = self.Mineral_list[self.inputnumber].Stix_Vp_Vs(self.P[i]*1e5,self.T[i])
                Vp.append(vp);Vs.append(vs);Rho.append(rho);K.append(k);G.append(g)
        if self.eos==2:
            for i in range(len(self.D)):
                vp,vs,rho,k,g = self.Mineral_list[self.inputnumber].HP_Vp_Vs(self.P[i]*1e5,self.T[i])
                Vp.append(vp);Vs.append(vs);Rho.append(rho);K.append(k);G.append(g)
        if self.eos==3:
            for i in range(len(self.D)):
                vp,vs,rho,k,g = self.Mineral_list[self.inputnumber].Thermal_Vp_Vs(self.P[i]*1e5,self.T[i],order=2)
                Vp.append(vp);Vs.append(vs);Rho.append(rho) ;K.append(k);G.append(g)   
        if self.eos==4:
            for i in range(len(self.D)):
                vp,vs,rho,k,g = self.Mineral_list[self.inputnumber].Thermal_Vp_Vs(self.P[i]*1e5,self.T[i],order=3)
                Vp.append(vp);Vs.append(vs);Rho.append(rho)  ;K.append(k);G.append(g)
        if self.eos==5:
            for i in range(len(self.D)):
                vp,vs,rho,k,g = self.Mineral_list[self.inputnumber].Thermal_Vp_Vs(self.P[i]*1e5,self.T[i],order=4)
                Vp.append(vp);Vs.append(vs);Rho.append(rho)  ;K.append(k);G.append(g)
                


        if self.plotoption ==2:  
            self.ax.cla()
            self.ax.plot(self.D,Vp,'g',label='Vp')
            self.ax.set_ylabel('$V_p$ (km/s)', color='k', fontname="Times New Roman",fontsize=10)
        
        if self.plotoption ==3:    
            self.ax.cla()
            self.ax.plot(self.D,Vs,'r',label='Vs')
            self.ax.set_ylabel('$V_s$ (km/s)', color='k', fontname="Times New Roman",fontsize=10)
            
        if self.plotoption ==4:  
            self.ax.cla()
            self.ax.plot(self.D,Vp,'g',label='Vp')
            self.ax.plot(self.D,Vs,'r',label='Vs')
            self.ax.set_ylabel('V (km/s)', color='k', fontname="Times New Roman",fontsize=10)
            
        if self.plotoption ==5:      
            self.ax.cla()
            self.ax.plot(self.D,Rho,'b',label='Rho')
            self.ax.set_ylabel('Density (kg/m$^3$)', color='k', fontname="Times New Roman",fontsize=10)

        if self.plotoption ==6:      
            self.ax.cla()
            self.ax.plot(self.D,np.array(K)/1e9,'b',label='K')
            self.ax.set_ylabel('Bulk Modulus (GPa)', color='k', fontname="Times New Roman",fontsize=10)
 
        if self.plotoption ==7:     
            self.ax.cla()
            self.ax.plot(self.D,np.array(G)/1e9,'b',label='G')
            self.ax.set_ylabel('Shear Modulus (GPa)', color='k', fontname="Times New Roman",fontsize=10)
            
        self.ax.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
        self.ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, borderaxespad=0.)
        self.ax.set_xlim(10,700)
        self.qmc.draw()       
        self.Vp=Vp
        self.Vs=Vs
        self.Rho=Rho
        self.K = K
        self.G = G

   
    
    def Plot_layout(self):
        self.Plot_Widget=QWidget()
        
        self.layout_fig = QVBoxLayout(self.Plot_Widget)
        self.figure = Figure(dpi=200)
        self.ax = self.figure.add_subplot(111) 
        self.qmc = FigureCanvas(self.figure)  
        #self.qmc.mpl_connect('button_press_event', self.OnClick)
        self.ntb = NavigationToolbar(self.qmc,self)      
        self.layout_fig.addWidget(self.qmc)
        self.layout_fig.addWidget(self.ntb)
        
        self.Export = QPushButton()
        self.Export.setText('Export')
        self.Export.setAutoDefault(False)
        self.Export.clicked.connect(self.export)
        self.layout_fig.addWidget(self.Export)
        
        self.Plot_Widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
    def export(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)        
        print (fileName)
        if fileName:
            f = open(fileName,'w')
            print ('open',fileName)
            f.write('Depth  Pressure  Temperature  Vp  Vs  Rho K G')
            f.write('\n')
            for i in range(len(self.P)):
                f.write(str(self.D[i]))
                f.write(' ' )        
                f.write(str(self.P[i]))
                f.write(' ' )
                f.write(str(self.T[i]))
                f.write(' ' )
                f.write(str(self.Vp[i]))   
                f.write(' ' )
                f.write(str(self.Vs[i])) 
                f.write(' ' )
                f.write(str(self.Rho[i]))                 
                f.write('\n')
                f.write(str(self.K[i]))                 
                f.write('\n')
                f.write(str(self.G[i]))                 
                f.write('\n')                
            f.close()
            print ('save')
            


if __name__ == "__main__":
    D=[];P=[];T=[]
    file = open('PT.txt')
    for line in file:
        line=line.split()
        D.append(float(line[0]))
        P.append(float(line[1]))
        T.append(float(line[2]))
    file.close()
    qapp = QApplication(sys.argv)
    GUI = GUI_user_input(D=D,P=P,T=T)
    #GUI = GUI_table()
    GUI.show()


#==============================================================================
#     Pressure= wfQ_ApplicationWindow(name="Pressure (Bar)")
#     Pressure.change(x=[1,2,3],y=[1,2,3])    
#     Pressure.exec_()
#==============================================================================
    sys.exit(qapp.exec_())    