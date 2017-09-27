# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 15:59:28 2016

@author: Fei

Update 10/3/2016
"""



from __future__ import print_function
import sys
import os
import re
import numpy as np
from scipy.odr import ODR, Model, RealData
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

try:
    from PyQt4.QtCore import (PYQT_VERSION_STR, QFile, QFileInfo, QSettings,
            QT_VERSION_STR, QTimer, QVariant, Qt, SIGNAL)
    from PyQt4.QtGui import (QAction, QActionGroup, QApplication,QSizePolicy,
            QDockWidget, QFileDialog, QFrame, QIcon, QImage, QImageReader,QColor,
            QImageWriter, QInputDialog, QKeySequence, QLabel, QListWidget,
            QMainWindow, QMessageBox, QPainter, QPixmap, QPrintDialog,QSplitter,QLineEdit,QComboBox,QCheckBox,
            QPrinter, QSpinBox,QStandardItemModel,QWidget,QGridLayout,QTableView,QDialog,QStandardItem,QPushButton)
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar    
except:
    from PyQt5.QtWidgets import (QSizePolicy,QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar   
    
try:   
    from .textEditor import TextEditor
except:
    from textEditor import TextEditor

# This define text formation
font = {'family' : 'serif',  
    'color'  : 'blue',  
    'weight' : 'normal',  
    'size'   : 10,  
    } 


class FigureCanvas3D(FigureCanvas):
    '''
    This class modify FigureCanvas from matplotlib package, and modify to show 3D plots
    '''
    def __init__(self,input1):
        X,Y,Z,a,b,c,d,strz,strx = input1
        self.fig = Figure(dpi=70)
        #self.fig.suptitle("this is  the  figure  title",  fontsize=12)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
            QSizePolicy.Expanding,
            QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.ax = Axes3D(self.fig) # Canvas figure must be created for mouse rotation
        self.ax.scatter3D(X,Y,Z,c='black', picker=5)           
        X1=np.arange(0.9*np.amin(X),1.1*np.amax(X),0.1)
        Y1=np.arange(0.9*np.amin(Y),1.1*np.amax(Y),0.1)
        X1,Y1= np.meshgrid(X1, Y1,sparse=True)
        Z1=b*X1+c*Y1+d*X1*Y1+a
        self.ax.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
        self.ax.plot_wireframe(X1,Y1,Z1,color='b')
        self.ax.set_xlabel("Water Content (wt%)",fontdict=font)
        self.ax.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
        self.ax.set_zlabel(strz,fontdict=font)
        self.fig.suptitle(strz,fontdict=font).set_fontsize(15)
        self.ax.set_xlim(0.9*np.amin(X),1.1*np.amax(X))
        self.ax.set_ylim(0.9*np.amin(Y),1.1*np.amax(Y))
        self.ax.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))
        self.ax.text2D(0.2,0.7,strx,transform=self.ax.transAxes)
        aa = np.linspace(0.9*np.amin(Z1), 1.1*np.amax(Z1), num=4)
        self.ax.set_zticks(aa)
        
        self.format_coord_org = self.ax.format_coord
        self.ax.format_coord = self.report_pixel
    
    def report_pixel(self, xd, yd):
        s = self.format_coord_org(xd, yd)
        s = s.replace(",", " ")
        return s
    
    def update_view(self):
        for pltt in self.plots:
            self.ax.collections.remove(pltt)
            self.plots = []
            self.draw()
            

 
    
 

class Regression(object):
    """
    This class do a linear regression
    """
    def __init__(self,name=None,water=None,Iron=None,flag=None):
        self.name=name
        self.key,self.content = self.read_data()
        self.Change=False
        if flag is not None:
            for i in range(len(flag)):
                if i> len(self.Flag):
                    break
                self.Flag[i]=float(flag[i])
        if self.name == "Olivine":
            self.pressure_deri = [4.217,1.462]
        if self.name == 'Wadsleyite':
            self.pressure_deri = [4.322,1.444]
        if self.name == 'Ringwoodite':
            self.pressure_deri = [4.22,1.354]

        try:
            self.water_content_error=self.content['H2O_err']
            self.iron_content_error=self.content['iron_err']
            self.K_error=self.content['K_err']
            self.G_error=self.content['G_err']
            self.Rho_error=self.content['Rho_err']
        except:
            self.water_content_error=np.zeros(self.number_of_data)+0.0
            self.iron_content_error=np.zeros(self.number_of_data)+0.0
            self.K_error=np.zeros(self.number_of_data)+0.0
            self.G_error=np.zeros(self.number_of_data)+0.0
            self.Rho_error=np.zeros(self.number_of_data) +0.0         
 
    def change_flag(self,flag):
        self.Flag = flag
        
    def original_flag(self):
        if self.name == "Olivine":
            self.Flag = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        if self.name == 'Wadsleyite':
            self.Flag = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        if self.name == 'Ringwoodite':
            self.Flag = [1,1,1,1,1,1,1,1,1,1,1,1,1]
        
    
    def Regression_Plane(self,x,x_error,y,y_error,z,z_error):
        """
        Do the regression, here use ODR package from Scipy which use orthogonal distance regression
        """
        data=RealData([x,y],z)#,sx=[x_error,y_error],sy=z_error)
        def func(beta,data):
            x,y = data
            a,b,c = beta
            return a*x+b*y+c    
        model = Model(func)
        odr = ODR(data, model,[100,100,100])
        res = odr.run()
        return res.beta, res.sd_beta    
    
    
    def Clean_content(self):
        #Clean all content
        self.Flag = []
        self.water_content= []
        self.iron_content= []
        self.K= []
        self.G= []
        self.Rho= []
        self.reference =  []
        
    def Return_original_data(self,num):
        #return data directly from txt file
        return [self.Flag[num],self.water_content[num],self.iron_content[num],self.K[num],self.G[num],self.Rho[num],self.reference[num]]

        
    def Change_data_from_table(self,flag):
        # user change data in table view
        # and here return data from table instead from txt file
        self.key,self.content = self.read_data()
        for i in range(len(flag)):
            if i> len(self.Flag):
                break
            self.Flag[i]=float(flag[i])
        for i in range(len(self.Flag)):
            if self.Flag[i] == 0:
                del self.water_content[i]
                del self.iron_content[i]
                del self.K[i]
                del self.G[i]
                del self.Rho[i]
                del self.reference[i]
                self.number_of_data = len(self.K)
        return None

    def Return(self):
        # based on user change return different data
        if self.Change == False:
            # return data without user input
            return self.function_K(),self.function_G(),self.function_Rho()
        else:
            # return data with user input 
            return self.UserK,self.UserG,self.UserRho

        
    def read_data(self):
        # read data from txt file
        try:
            address=os.path.join(os.path.dirname(__file__),'EXPDATA',self.name+'.txt')
            self.address=address
            file=open(address,'r+')
        except:
            address=os.path.join(os.path.dirname(__file__),'Mineral_Physics','EXPDATA',self.name+'.txt')
            self.address=address
            file=open(address,'r+')           
        for i in file:
            name=i.split(',')
            break
        dictionary=dict()  
        for i in name:   
            dictionary[i]=[]    
        for line in file:
           line=line.split(',')
           #print (line)
           for i in range(len(line)):  
               try:
                   dictionary[name[i]].append(float(line[i]))
               except:
                   dictionary[name[i]].append((line[i])) 
        file.close()
        self.Flag = dictionary['Flag']
        self.water_content=dictionary['H2O'];
        self.iron_content=dictionary['Iron']
        self.K=dictionary['K']
        self.G=dictionary['G']
        self.Rho=dictionary['Rho']
        self.reference = dictionary['Author']
        self.number_of_data=len(self.K)        
        return name,dictionary  
    


    def func(self,beta,data):
        # linear function use to do the regression 
        x,y = data
        a,b,c = beta
        return a*x+b*y+c

    '''
    3 functions that return regression coefficient 
    '''
    def function_K(self):
        return self.Regression_Plane(self.water_content,self.water_content_error,
                    self.iron_content,self.iron_content_error,self.K,self.K_error)     
    def function_G(self):
        return self.Regression_Plane(self.water_content,self.water_content_error,
                    self.iron_content,self.iron_content_error,self.G,self.G_error) 
    def function_Rho(self):
        return self.Regression_Plane(self.water_content,self.water_content_error,
                    self.iron_content,self.iron_content_error,self.Rho,self.Rho_error) 

    
    def Show_fit_function(self,beta,sd_beta,name,error=True):
        # function that show fit results.             
        b11=beta[0];c11=beta[1];a11=beta[2];
        b12=abs(sd_beta[0]);c12=abs(sd_beta[1]);a12=abs(sd_beta[2]);
        if error == True:
            string1=name + ' %5.2f(%5.2f)' %(a11,a12)
            def plus_Minus(b11,b12):
                if b11>=0:
                    return '+%5.2f(%5.2f)' %(b11,b12)
                else:
                    return '%5.2f(%5.2f)' %(b11,b12)
            string2=plus_Minus(b11,b12)
            string3=plus_Minus(c11,c12)
            
            string=string1+string2+'*CWater'+string3+'*XFe'+'\n'
            
        else:
            string1='%.2f' %(a11)
            def plus_Minus(b11,b12):
                if b11>=0:
                    return '+%.2f' %(b11)
                else:
                    return '%.2f' %(b11)
            string2=plus_Minus(b11,b12)
            string3=plus_Minus(c11,c12)
            if name == "K'":
                string=string2+'*CWater'+string3+'*XFe+'+string1+","+name+"="+str(self.pressure_deri[0])     
            elif name == "G'":
                string=string2+'*CWater'+string3+'*XFe+'+string1+","+name+"="+str(self.pressure_deri[1])                   
            else:
                string=string2+'*CWater'+string3+'*XFe+'+string1+""+name+""+''                     
        return string 
        
    
        
    def PLOT(self,return_fig=True):  
        '''
        This function return PLOT using data from table 
        '''
        self.Change = False
        if self.name == "Olivine":
            self.pressure_deri = [4.217,1.462]
        if self.name == 'Wadsleyite':
            self.pressure_deri = [4.322,1.444]
        if self.name == 'Ringwoodite':
            self.pressure_deri = [4.22,1.354]
        self.str1=self.Show_fit_function(self.function_K()[0],self.function_K()[1],'')
        self.str2=self.Show_fit_function(self.function_G()[0],self.function_G()[1],'')
        self.str3=self.Show_fit_function(self.function_Rho()[0],self.function_Rho()[1],'')
                
        self.d1,self.b1,self.c1,self.a1=0,self.function_K()[0][0],self.function_K()[0][1],self.function_K()[0][2]
        self.d2,self.b2,self.c2,self.a2=0,self.function_G()[0][0],self.function_G()[0][1],self.function_G()[0][2]
        self.d3,self.b3,self.c3,self.a3=0,self.function_Rho()[0][0],self.function_Rho()[0][1],self.function_Rho()[0][2]
        
        if return_fig ==True:        
            fig1=plot(self.water_content,self.iron_content,self.K,self.a1,self.b1,self.c1,self.d1,self.name+" Bulk modulus",self.str1)
            fig2=plot(self.water_content,self.iron_content,self.G,self.a2,self.b2,self.c2,self.d2,self.name+" Shear modulus",self.str2)
            fig3=plot(self.water_content,self.iron_content,self.Rho,self.a3,self.b3,self.c3,self.d3,self.name+" Density",self.str3)
            return fig1,fig2,fig3
        else:    
            return [self.water_content,self.iron_content,self.K,self.a1,self.b1,self.c1,self.d1,self.name+" Bulk modulus",self.str1],[self.water_content,self.iron_content,self.G,self.a2,self.b2,self.c2,self.d2,self.name+" Shear modulus",self.str2],[self.water_content,self.iron_content,self.Rho,self.a3,self.b3,self.c3,self.d3,self.name+" Density",self.str3]
     
     
    def User_Change(self,K,G,Rho):
        self.Change = True
        self.UserK=K
        self.UserG = G
        self.UserRho = Rho
        self.pressure_deri = [float(self.UserK[0][-1]),float(self.UserG[0][-1])]
        #print ('mineral user change')
        

        
     
     
    def PLOT_input_formula(self,return_fig=True):  
        '''
        This function return PLOT using user input 
        '''
   
        self.str1=self.UserK[1]
        self.str2=self.UserG[1]
        self.str3=self.UserRho[1]
                
        self.b1,self.c1,self.a1=self.UserK[0][:3].astype(float)
        self.b2,self.c2,self.a2=self.UserG[0][:3].astype(float)
        self.b3,self.c3,self.a3=self.UserRho[0][:3].astype(float)
        self.d1,self.d2,self.d3 = 0,0,0
        
        #print ('pass parameters')
        
        if return_fig ==True:        
            fig1=plot(self.water_content,self.iron_content,self.K,self.a1,self.b1,self.c1,0,self.name+" Bulk modulus",self.str1)
            fig2=plot(self.water_content,self.iron_content,self.G,self.a2,self.b2,self.c2,0,self.name+" Shear modulus",self.str2)
            fig3=plot(self.water_content,self.iron_content,self.Rho,self.a3,self.b3,self.c3,0,self.name+" Density",self.str3)
            return fig1,fig2,fig3
        else:    
            return [self.water_content,self.iron_content,self.K,self.a1,self.b1,self.c1,self.d1,self.name+" Bulk modulus",self.str1],[self.water_content,self.iron_content,self.G,self.a2,self.b2,self.c2,self.d2,self.name+" Shear modulus",self.str2],[self.water_content,self.iron_content,self.Rho,self.a3,self.b3,self.c3,self.d3,self.name+" Density",self.str3]
         
             
     
         
class Regression_PLOT_PyQt(QDialog):
    """
    This calss return a PyQt GUI with a regression plot
    """       
    def __init__(self,Minerals=None,string=None,flag = None):
        super(Regression_PLOT_PyQt, self).__init__()
        if string is not None:
            self.stringK,self.stringG,self.stringRho,self.user_input = string
        else:
            self.stringK=None;self.stringG=None;self.stringRho=None;self.user_input=False
        self.resize(1400, 600)
        self.Minerals = Minerals

        self.dirty=False
        #self.user_input = user_input

        self.Minerals.original_flag()
        self.Minerals.read_data()
        self.table = QTableView()
        self.model = QStandardItemModel(25,7,self)
        self.model.setHorizontalHeaderLabels(['flag','Water Content','Iron Content','K (Gpa)','G (Gpa)', 'Rho (g/cm³)','Reference'])
        for i in range(len(self.Minerals.Flag)):
                a=self.Minerals.Return_original_data(i)
                for j in range(0,7):
                    item = QStandardItem(str(a[j]))
                    self.model.setItem(i, j,item)  
                    if j != 0:
                        item.setFlags(Qt.ItemIsEnabled)
                        item.setBackground(QColor(211,211,211))

        if flag is not None:
            self.Minerals.change_flag(flag)
            for i in range(len(self.Minerals.Flag)):
                item = QStandardItem(str(self.Minerals.Flag[i]))
                self.model.setItem(i, 0,item)                 
            
        self.table.setModel(self.model)
        
        self.button = QPushButton('Update and use in thermoelastic model')
        self.button.clicked.connect(self.Update)  
        self.button.setAutoDefault(False)
        self.button1 = QPushButton('Add data file ')
        self.button1.clicked.connect(self.Export)   
        self.button1.setAutoDefault(False)        
        self.layout = QGridLayout()


        
        self.label = QLabel()
        self.label.setText(        '''
        Please input equation, Water: water content (wt%) and Fe: iron content (mol%)
        for example -2.41*Water-30*Fe+81,K'=4.1,  
        ''')
        
        self.Kinput_formula = QLineEdit()
        self.Ginput_formula = QLineEdit()
        self.Rhoinput_formula = QLineEdit()
        if self.stringK is not None:
            self.Kinput_formula.setText(self.stringK)
            self.Ginput_formula.setText(self.stringG)
            self.Rhoinput_formula.setText(self.stringRho)
            self.Userinput()
        else:            
            self.Kinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_K()[0],self.Minerals.function_K()[1],"K'",error=False))
            self.Ginput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_G()[0],self.Minerals.function_G()[1],"G'",error=False))
            self.Rhoinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_Rho()[0],self.Minerals.function_Rho()[1],'',error=False))
        self.Kinput_formula.returnPressed.connect(self.Kformula)
        self.Ginput_formula.returnPressed.connect(self.Gformula)
        self.Rhoinput_formula.returnPressed.connect(self.Rhoformula)
        #self.connect(self.Kinput_formula,SIGNAL("returnPressed()"),self.Kformula)
        #self.connect(self.Ginput_formula,SIGNAL("returnPressed()"),self.Gformula)
        #self.connect(self.Rhoinput_formula,SIGNAL("returnPressed()"),self.Rhoformula)
        
        self.user_change_confrim = QPushButton('Change to user input')
        self.user_change_confrim.clicked.connect(self.Userinput)   
        self.user_change_confrim.setAutoDefault(False)

        
        self.extension = QWidget()
        self.extensionLayout = QGridLayout()
        self.extensionLayout.setContentsMargins(0, 0, 0, 0)

        labelK = QLabel()
        labelK.setText('K<sub>0</sub>=')
        labelG = QLabel()
        labelG.setText('G<sub>0</sub>=')
        labelRho = QLabel()
        labelRho.setText('Rho<sub>0</sub>=')

        self.extensionLayout.addWidget(labelK,0,0,1,3)
        self.extensionLayout.addWidget(labelG,1,0,1,3)
        self.extensionLayout.addWidget(labelRho,2,0,1,3)        
        self.extensionLayout.addWidget(self.Kinput_formula,0,1,1,3)
        self.extensionLayout.addWidget(self.Ginput_formula,1,1,1,3)
        self.extensionLayout.addWidget(self.Rhoinput_formula,2,1,1,3)
        self.extensionLayout.addWidget(self.user_change_confrim,3,0,9,4)
        
        
        self.extension.setLayout(self.extensionLayout)
       
        self.check_change = QCheckBox("user input")
        self.check_change.setChecked(False)
        self.check_change.toggled.connect(self.extension.setVisible)
        
        #self.PLOT(switch=True)
        self.PLOT()   
        self.layout.addWidget(self.table,0,1,1,17)
        self.layout.addWidget(self.button,2,0,1,9)
        self.layout.addWidget(self.button1,2,9,1,9)
        self.layout.addWidget(self.check_change,3,0,1,1)
        self.layout.addWidget(self.label,3,3,1,5)
        self.layout.addWidget(self.extension,4,0,1,9)        
        self.setLayout(self.layout)  
        
        self.extension.hide()

    def Kformula(self):
        self.stringK=str(self.Kinput_formula.text())
        self.K0 = np.array(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", self.stringK))
        return self.K0,self.stringK
        
    def Gformula(self):
        self.stringG=str(self.Ginput_formula.text())
        self.G0= np.array(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", self.stringG))
        return self.G0,self.stringG
        
    def Rhoformula(self):
        self.stringRho=str(self.Rhoinput_formula.text())
        self.Rho0 = np.array(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", self.stringRho))
        return self.Rho0,self.stringRho
    
    def Userinput(self):
        self.Minerals.User_Change(self.Kformula(),self.Gformula(),self.Rhoformula())
        self.user_input=True
        self.PLOT()
  

    
        
    def PLOT(self,switch=True):
        if self.dirty==False:
            #self.Minerals.Clean_content()
            params=[]
            for i in range(20):            
                index=self.model.index(i,0)
                try:
                    aa=(float(self.model.itemData(index)[0]))
                    #print (aa)
                    params.append(aa)
                except:
                    break
            #print (params)
            self.Minerals.Change_data_from_table(params)
            #print ('wf')
        else:          
            self.model.setHorizontalHeaderLabels(['flag','Water Content','Iron Content','K (Gpa)','G (Gpa)', 'Rho (g/cm3)','Reference'])
            for i in range(len(self.Minerals.Flag)):
                #a=self.Minerals.Return_original_data(i)
                item = QStandardItem(str(self.Minerals.Flag[i]))
                self.model.setItem(i, 0,item)  
#==============================================================================
#                 for j in range(1,7):
#                     item = QStandardItem(str(a[j]))
#                     self.model.setItem(i, j,item)  
#                     if j != 0:
#                         item.setFlags(Qt.ItemIsEnabled)
#                         item.setBackground(QColor(211,211,211))
#==============================================================================
                        
        #print (self.Minerals.Change)     
            
            
        if self.user_input == False:
            self.a,self.b,self.c = self.Minerals.PLOT(return_fig=False)       
            #print (self.Minerals.number_of_data)
            self.Kinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_K()[0],self.Minerals.function_K()[1],"K'",error=False))
            self.Ginput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_G()[0],self.Minerals.function_G()[1],"G'",error=False))
            self.Rhoinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_Rho()[0],self.Minerals.function_Rho()[1],'',error=False))

        else:
            self.a,self.b,self.c = self.Minerals.PLOT_input_formula(return_fig=False)   
       
            
            
        self.canvas1 = FigureCanvas3D(self.a)
        self.canvas2 = FigureCanvas3D(self.b)
        self.canvas3 = FigureCanvas3D(self.c)
        self.canvas1.mpl_connect('pick_event', self.onpick)
        self.canvas2.mpl_connect('pick_event', self.onpick)
        self.canvas3.mpl_connect('pick_event', self.onpick)
        
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.toolbar3 = NavigationToolbar(self.canvas3, self) 
 
        self.layout1_widget = QWidget()
        self.layout1 = QGridLayout(self.layout1_widget)
        self.layout1_widget.setFixedSize(600, 600)
        self.layout1.addWidget(self.canvas1,0,1,5,5)
        self.layout1.addWidget(self.toolbar1,5,1,1,5)       
        self.layout1.addWidget(self.canvas2,6,1,5,5)
        self.layout1.addWidget(self.toolbar2,11,1,1,5)  
        self.layout1.addWidget(self.canvas3,12,1,5,5)
        self.layout1.addWidget(self.toolbar3,17,1,1,5)       
        self.layout.addWidget(self.layout1_widget,0,0,1,1)
        
    def onpick(self,event):
        try:
            for i in range(6):
                self.model.item(self.ind,i+1).setBackground(QColor(211,211,211))
        except:
            pass
        
        count=-1
        for j in range(len((self.Minerals.Flag))):
            if self.Minerals.Flag[j]==1:
                count+=1
                if count == event.ind[0]:
                    self.ind=j
                    break
        #print (self.ind)       
        #self.ind = event.ind 
        for i in range(6):
            self.model.item(self.ind,i+1).setBackground(QColor(111,111,111))

        
    def Update(self):
        self.user_input = False
        self.dirty=False
        self.check_change.setChecked(False)
        self.PLOT()
    
    def Export(self):
        self.dirty=True
        dialog = TextEditor(name=self.Minerals.name)
        if dialog.exec_():
                pass
        self.Minerals.read_data()        
        self.PLOT()
        
    def ReturnString(self):
        aa = self.Minerals.Flag
        bb=str(int(aa[0]))
        for i in range(1,len(aa)):
            bb+=str(int(aa[i]))
        return [self.stringK,self.stringG,self.stringRho,self.user_input],bb



olivine=Regression('Olivine')
wadsleyte=Regression('Wadsleyite')
ringwoodite=Regression('Ringwoodite')        
        
        
if __name__ == "__main__":

#==============================================================================
#     qApp = QApplication(sys.argv)
#     app = Regression_PLOT_PyQt(olivine,["-4.23*CWater+9.96*XFe+128.27,K'=4.217",
#  "-2.36*CWater-30.20*XFe+81.14,G'=4.217",
#  '-301.14*CWater+142.69*XFe+3313.46',
#  True],flag='1111111111111111')
#     app.show()
#     sys.exit(qApp.exec_())        
#==============================================================================
    font = {'family' : 'serif',  
        'color'  : 'blue',  
        'weight' : 'normal',  
        'size'   : 10,  
        } 
    Minerals = olivine
    Minerals.read_data()
    a,b,c = Minerals.PLOT(return_fig=False)       
    

    Xa,Ya,Za,aa,ba,ca,da,strza,strxa = a
    Xb,Yb,Zb,ab,bb,cb,db,strzb,strxb = b
    Xc,Yc,Zc,ac,bc,cc,dc,strzc,strxc = c    
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(331, projection='3d')
    ax2 = fig.add_subplot(332, projection='3d')
    ax3 = fig.add_subplot(333, projection='3d')
    ax4 = fig.add_subplot(334, projection='3d')
    ax5 = fig.add_subplot(335, projection='3d')
    ax6 = fig.add_subplot(336, projection='3d')    
    ax7 = fig.add_subplot(337, projection='3d')
    ax8 = fig.add_subplot(338, projection='3d')
    ax9 = fig.add_subplot(339, projection='3d')    
    
    
    ax1.scatter3D(Xa,Ya,Za,c='black', picker=5)      
    ax2.scatter3D(Xb,Yb,Zb,c='black', picker=5)      
    ax3.scatter3D(Xc,Yc,Zc,c='black', picker=5) 
     
    X1=np.arange(0.9*np.amin(Xa),1.1*np.amax(Xa),0.1)
    Y1=np.arange(0.9*np.amin(Ya),1.1*np.amax(Ya),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=ba*X1+ca*Y1+da*X1*Y1+aa
    ax1.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax1.plot_wireframe(X1,Y1,Z1,color='b')
    ax1.set_xlabel("Water Content (wt%)",fontdict=font)
    ax1.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax1.set_zlabel(strza,fontdict=font)
    #fig.suptitle(strza,fontdict=font).set_fontsize(15)
    ax1.set_xlim(0.9*np.amin(Xa),1.1*np.amax(Xa))
    ax1.set_ylim(0.9*np.amin(Ya),1.1*np.amax(Ya))
    ax1.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))
    ax1.set_zlim3d(120, 140)
    ax1.set_zticks([120, 130,140])
    
    #ax1.text2D(0,0,strxa)
    X1=np.arange(0.9*np.amin(Xb),1.1*np.amax(Xb),0.1)
    Y1=np.arange(0.9*np.amin(Yb),1.1*np.amax(Yb),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=bb*X1+cb*Y1+db*X1*Y1+ab
    ax2.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax2.plot_wireframe(X1,Y1,Z1,color='b')
    ax2.set_xlabel("Water Content (wt%)",fontdict=font)
    ax2.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax2.set_zlabel(strzb,fontdict=font)
    #fig.suptitle(strzb,fontdict=font).set_fontsize(15)
    ax2.set_xlim(0.9*np.amin(Xb),1.1*np.amax(Xb))
    ax2.set_ylim(0.9*np.amin(Yb),1.1*np.amax(Yb))
    ax2.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))    
    ax2.set_zticks([50,60,70,80])


    
    #ax1.text2D(0,0,strxa)
    
    X1=np.arange(0.9*np.amin(Xc),1.1*np.amax(Xc),0.1)
    Y1=np.arange(0.9*np.amin(Yc),1.1*np.amax(Yc),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=bc*X1+cc*Y1+dc*X1*Y1+ac
    ax3.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax3.plot_wireframe(X1,Y1,Z1,color='b')
    ax3.set_xlabel("Water Content (wt%)",fontdict=font)
    ax3.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax3.set_zlabel(strzc,fontdict=font)
    #fig.suptitle(strzc,fontdict=font).set_fontsize(15)
    ax3.set_xlim(0.9*np.amin(Xc),1.1*np.amax(Xc))
    ax3.set_ylim(0.9*np.amin(Yc),1.1*np.amax(Yc))
    ax3.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))     
    ax3.set_zticks([2800,3200,3600,4000])
    
    
    

    '''
    WAD
    '''
    Minerals = wadsleyte
    Minerals.read_data()
    a,b,c = Minerals.PLOT(return_fig=False)       
    

    Xa,Ya,Za,aa,ba,ca,da,strza,strxa = a
    Xb,Yb,Zb,ab,bb,cb,db,strzb,strxb = b
    Xc,Yc,Zc,ac,bc,cc,dc,strzc,strxc = c    
    
    #f, (ax1, ax2, ax3) = plt.subplots(3,1, projection='3d')
    
    ax4.scatter3D(Xa,Ya,Za,c='black', picker=5)      
    ax5.scatter3D(Xb,Yb,Zb,c='black', picker=5)      
    ax6.scatter3D(Xc,Yc,Zc,c='black', picker=5) 
     
    X1=np.arange(0.9*np.amin(Xa),1.1*np.amax(Xa),0.1)
    Y1=np.arange(0.9*np.amin(Ya),1.1*np.amax(Ya),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=ba*X1+ca*Y1+da*X1*Y1+aa
    ax4.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax4.plot_wireframe(X1,Y1,Z1,color='b')
    ax4.set_xlabel("Water Content (wt%)",fontdict=font)
    ax4.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax4.set_zlabel(strza,fontdict=font)
    #fig.suptitle(strza,fontdict=font).set_fontsize(15)
    ax4.set_xlim(0.9*np.amin(Xa),1.1*np.amax(Xa))
    ax4.set_ylim(0.9*np.amin(Ya),1.1*np.amax(Ya))
    ax4.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))
    ax4.set_zlim3d(120, 200)
    ax4.set_zticks([120,140,160,180,200])
    
    #ax1.text2D(0,0,strxa)
    X1=np.arange(0.9*np.amin(Xb),1.1*np.amax(Xb),0.1)
    Y1=np.arange(0.9*np.amin(Yb),1.1*np.amax(Yb),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=bb*X1+cb*Y1+db*X1*Y1+ab
    ax5.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax5.plot_wireframe(X1,Y1,Z1,color='b')
    ax5.set_xlabel("Water Content (wt%)",fontdict=font)
    ax5.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax5.set_zlabel(strzb,fontdict=font)
    #fig.suptitle(strzb,fontdict=font).set_fontsize(15)
    ax5.set_xlim(0.9*np.amin(Xb),1.1*np.amax(Xb))
    ax5.set_ylim(0.9*np.amin(Yb),1.1*np.amax(Yb))
    ax5.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))    
    ax5.set_zticks([80,100,120])


    
    #ax1.text2D(0,0,strxa)
    
    X1=np.arange(0.9*np.amin(Xc),1.1*np.amax(Xc),0.1)
    Y1=np.arange(0.9*np.amin(Yc),1.1*np.amax(Yc),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=bc*X1+cc*Y1+dc*X1*Y1+ac
    ax6.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax6.plot_wireframe(X1,Y1,Z1,color='b')
    ax6.set_xlabel("Water Content (wt%)",fontdict=font)
    ax6.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax6.set_zlabel(strzc,fontdict=font)
    #fig.suptitle(strzc,fontdict=font).set_fontsize(15)
    ax6.set_xlim(0.9*np.amin(Xc),1.1*np.amax(Xc))
    ax6.set_ylim(0.9*np.amin(Yc),1.1*np.amax(Yc))
    ax6.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))     
    ax6.set_zticks([3000,3300,3600,3900])
    
    
    
#==============================================================================
#     fig1.savefig('Test4.jpg',dpi=600)
#     fig2.savefig('Test5.jpg',dpi=600)
#     fig3.savefig('Test6.jpg',dpi=600)
#==============================================================================

    
    '''
    RI
    '''
    Minerals = ringwoodite
    Minerals.read_data()
    a,b,c = Minerals.PLOT(return_fig=False)       
    

    Xa,Ya,Za,aa,ba,ca,da,strza,strxa = a
    Xb,Yb,Zb,ab,bb,cb,db,strzb,strxb = b
    Xc,Yc,Zc,ac,bc,cc,dc,strzc,strxc = c    
    

    #f, (ax1, ax2, ax3) = plt.subplots(3,1, projection='3d')
    
    ax7.scatter3D(Xa,Ya,Za,c='black', picker=5)      
    ax8.scatter3D(Xb,Yb,Zb,c='black', picker=5)      
    ax9.scatter3D(Xc,Yc,Zc,c='black', picker=5) 
     
    X1=np.arange(0.9*np.amin(Xa),1.1*np.amax(Xa),0.1)
    Y1=np.arange(0.9*np.amin(Ya),1.1*np.amax(Ya),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=ba*X1+ca*Y1+da*X1*Y1+aa
    ax7.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax7.plot_wireframe(X1,Y1,Z1,color='b')
    ax7.set_xlabel("Water Content (wt%)",fontdict=font)
    ax7.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax7.set_zlabel(strza,fontdict=font)
    #fig.suptitle(strza,fontdict=font).set_fontsize(15)
    ax7.set_xlim(0.9*np.amin(Xa),1.1*np.amax(Xa))
    ax7.set_ylim(0.9*np.amin(Ya),1.1*np.amax(Ya))
    ax7.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))
    ax7.set_zticks([140,165,190,215])
    
    #ax1.text2D(0,0,strxa)
    X1=np.arange(0.9*np.amin(Xb),1.1*np.amax(Xb),0.1)
    Y1=np.arange(0.9*np.amin(Yb),1.1*np.amax(Yb),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=bb*X1+cb*Y1+db*X1*Y1+ab
    ax8.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax8.plot_wireframe(X1,Y1,Z1,color='b')
    ax8.set_xlabel("Water Content (wt%)",fontdict=font)
    ax8.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax8.set_zlabel(strzb,fontdict=font)
    #fig.suptitle(strzb,fontdict=font).set_fontsize(15)
    ax8.set_xlim(0.9*np.amin(Xb),1.1*np.amax(Xb))
    ax8.set_ylim(0.9*np.amin(Yb),1.1*np.amax(Yb))
    ax8.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))    
    ax8.set_zticks([45,70,95,120])


    
    #ax1.text2D(0,0,strxa)
    
    X1=np.arange(0.9*np.amin(Xc),1.1*np.amax(Xc),0.1)
    Y1=np.arange(0.9*np.amin(Yc),1.1*np.amax(Yc),0.1)
    X1,Y1= np.meshgrid(X1, Y1,sparse=True)
    Z1=bc*X1+cc*Y1+dc*X1*Y1+ac
    ax9.plot_surface(X1,Y1,Z1,color='b',alpha=0.1)
    ax9.plot_wireframe(X1,Y1,Z1,color='b')
    ax9.set_xlabel("Water Content (wt%)",fontdict=font)
    ax9.set_ylabel("Iron Content (mol%)",fontdict=font)#.set_fontsize(10)
    ax9.set_zlabel(strzc,fontdict=font)
    #fig.suptitle(strzc,fontdict=font).set_fontsize(15)
    ax9.set_xlim(0.9*np.amin(Xc),1.1*np.amax(Xc))
    ax9.set_ylim(0.9*np.amin(Yc),1.1*np.amax(Yc))
    ax9.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))     
    ax9.set_zticks([3000,4000,5000])
    
    
    fig.savefig('1.eps',dpi=200)
    
    
    
    

        
        
        
        
        
        
        
        
        


    

