# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 15:59:28 2016
This class contian GUI to do the regression and retun results
@author: Fei
Update 10/3/2016
"""



from __future__ import print_function
import sys
import os
import re
import numpy as np
from scipy.odr import ODR, Model, RealData,Data
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import functools
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
    def __init__(self,input1,error,error_show=True):
        x,y,z,a,b,c,d,strz,strx = input1
        
        x_error,y_error,z_error = error
        
        X=[];Y=[];Z=[];X_error=[];Y_error=[];Z_error=[];
        for i in range(len(z)):
            try:
                Z.append(float(z[i]))
                Z_error.append(float(z_error[i]))
                X.append(float(x[i]))
                X_error.append(float(x_error[i]))
                Y.append(float(y[i]))
                Y_error.append(float(y_error[i]))
            except:
                pass                
#==============================================================================
#             try:
#                 z[i] = float(z[i])
#             except:
#                 pass
#             if isinstance(z[i], float):
#                 Z.append(float(z[i]))
#                 X.append(float(x[i]))
#                 Y.append(float(y[i]))
#             else:
#                 pass
#==============================================================================
                
        lingspace = 2
        self.fig = Figure(dpi=70)
        #self.fig.suptitle("this is  the  figure  title",  fontsize=12)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
            QSizePolicy.Expanding,
            QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.ax = Axes3D(self.fig) # Canvas figure must be created for mouse rotation
        #self.ax =self.fig.gca(projection='3d')
        #self.ellipsoid()
        if error_show == True:
            for i in range(len(X)):
                self.ellipsoid(a=X[i],b=Y[i],c=Z[i],a_err=X_error[i],b_err=Y_error[i],c_err=Z_error[i])
        else:
            pass
            
        self.ax.scatter3D(X,Y,Z,c='b', picker=5)   
        #print (X,Y,Z)        
        X1=np.arange(0.9*np.amin(X),1.1*np.amax(X),0.1)
        Y1=np.arange(0.9*np.amin(Y),1.2*np.amax(Y),0.05)
        X1,Y1= np.meshgrid(X1, Y1,sparse=True)
        Z1=b*X1+c*Y1+d*X1*Y1+a
        self.ax.plot_surface(X1,Y1,Z1,color='r',alpha=0.1)
        #self.ax.plot_wireframe(X1,Y1,Z1,color='b')
        self.ax.set_xlabel("\nWater Content (wt%)",fontdict=font,linespacing=lingspace)
        self.ax.set_ylabel("\nIron Content (mol%)",fontdict=font,linespacing=lingspace)#.set_fontsize(10)
        self.ax.set_zlabel(strz,fontdict=font,linespacing=lingspace)
        self.fig.suptitle(strz,fontdict=font).set_fontsize(15)
        self.ax.set_xlim(0.9*np.amin(X),1.1*np.amax(X))
        self.ax.set_ylim(0.9*np.amin(Y),1.1*np.amax(Y))
        self.ax.set_zlim(0.9*np.amin(Z1),1.1*np.amax(Z1))
        self.ax.text2D(0.2,0.7,strx,transform=self.ax.transAxes)
        #self.ax.xaxis._axinfo['label']['space_factor'] = 20
        
        
        #self.format_coord_org = self.ax.format_coord
        #self.ax.format_coord = self.report_pixel
        
    
    def report_pixel(self, xd, yd):
        s = self.format_coord_org(xd, yd)
        s = s.replace(",", " ")
        return s
    
    def update_view(self):
        for pltt in self.plots:
            self.ax.collections.remove(pltt)
            self.plots = []
            self.draw()
            
    
    def ellipsoid(self,a=1,b=1,c=1,a_err=0.01,b_err=0.0001,c_err=1):

        
        u = np.linspace(0, 2 * np.pi, 10)
        v = np.linspace(0, np.pi, 10)
        
#==============================================================================
#         if b_err <=0.001:
#             b_err=0.01
#==============================================================================
        # Cartesian coordinates that correspond to the spherical angles:
        # (this is the equation of an ellipsoid):
        x = a_err * np.outer(np.cos(u), np.sin(v))+a
        y = b_err * np.outer(np.sin(u), np.sin(v))+b
        z = c_err * np.outer(np.ones_like(u), np.cos(v))+c
        
        # Plot:
        self.ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')
 
    
 

class Regression(object):
    """
    This class do a linear regression
    """
    def __init__(self,name=None,water=None,Iron=None,flag=None):
        self.name=name
        self.key,self.content = self.read_data()
        self.Change=False
        self.regression_methods = 1
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
        
    
    def Regression_Plane(self,x,x_error,y,y_error,z,z_error,methods=2):
        """
        Do the regression, here use ODR package from Scipy which use orthogonal distance regression
        """
        X=[];Y=[];Z=[];X_error=[];Y_error=[];Z_error=[]
        for i in range(len(z)):
            try:
                Z.append(float(z[i]))
                X.append(float(x[i]))
                Y.append(float(y[i]))
                X_error.append(float(x_error[i]))
                Y_error.append(float(y_error[i]))
                Z_error.append(float(z_error[i]))
            except:
                pass 
            
        if methods==2:
            #data=RealData([X,Y],Z)#,sx=[x_error,y_error],sy=z_error)
            data = Data([X,Y],Z,wd=[1,1])
            def func(beta,data):
                x,y = data
                a,b,c = beta
                return a*x+b*y+c    
            model = Model(func)
            odr = ODR(data, model,[100,100,100])
            res = odr.run()
            #print ('dond odr')
            return res.beta, res.sd_beta    
        
        elif methods ==0 :
            
            def plane_method(x2, y2, params):
                a = params[0]
                b = params[1]
                c = params[2]
                z = a*x2 + b*y2 + c
                #print (x2,y2)
                return z
        
            def error_method(params, points):
                result = []
                for (x1,y1,z1) in points:
                    #up = np.abs(params[0]*x+params[1]*y -z +params[2])
                    #down = np.sqrt(params[0]**2 +params[1]**2 +1)
                    plane_z = plane_method(x1, y1, params)
                    #print (plane_z)
                    diff = abs(plane_z - z1)
                    #diff = up/down
                    result .append( diff**2)
                return result

            points=[]
            for i in range(len(Z)):
                points.append([X[i],Y[i],Z[i]])
            #print (points)

            fun = functools.partial(error_method, points=points)
            params0 = [100, 100, 100]
            res = leastsq(fun, params0,full_output=1)
            cov = res[1]
            perr = np.diag(cov)
            #print (perr)
            a = res[0][0]
            b = res[0][1]
            c = res[0][2]   
            return np.array([a,b,c]),np.array(perr)
        
        elif methods ==1:

            data=RealData([X,Y],Z,sx=[X_error,Y_error],sy=Z_error)
            #data = Data([X,Y],Z,wd=[1,100])
            def func(beta,data):
                x,y = data
                a,b,c = beta
                return a*x+b*y+c    
            model = Model(func)
            odr = ODR(data, model,[100,100,100])
            res = odr.run()
           # print ('dond odr eror')
            return res.beta, res.sd_beta  
        
        elif methods==3:
            data = Data([X,Y],Z,wd=[1,100])
            def func(beta,data):
                x,y = data
                a,b,c = beta
                return a*x+b*y+c    
            model = Model(func)
            odr = ODR(data, model,[100,100,100])
            res = odr.run()
            #print ('dond odr')
            return res.beta, res.sd_beta          
        
        else:
            raise ('wrong regression')
            
        
    
    
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
        return [self.Flag[num],self.dictionary['H2O'][num],self.dictionary['Iron'][num],self.dictionary['K'][num],self.dictionary['G'][num],self.dictionary['Rho'][num],self.reference[num]]

        
    def Change_data_from_table(self,flag):
        # user change data in table view
        # and here return data from table instead from txt file
        self.key,self.content = self.read_data()
        for i in range(len(flag)):
            if i> len(self.Flag):
                break
            self.Flag[i]=float(flag[i])
        del_num = 0
        
        self.water_content = list(self.water_content)
        self.iron_content= list(self.iron_content)
        self.K = list(self.K)
        self.G = list(self.G)
        self.Rho = list(self.Rho)
        self.reference = list(self.reference)
        
        for i in range(len(self.Flag)):
            if self.Flag[i] == 0:
                del self.water_content[i-del_num]
                del self.iron_content[i-del_num]
                del self.K[i-del_num]
                del self.G[i-del_num]
                del self.Rho[i-del_num]
                del self.reference[i-del_num]
                del_num  += 1
        
        self.number_of_data = len(self.K)
        return None

    def Return(self,methods=2):
        # based on user change return different data
        if self.Change == False:
            # return data without user input
            return self.function_K(methods = self.regression_methods),self.function_G(methods = self.regression_methods),self.function_Rho(methods = self.regression_methods)
        else:
            # return data with user input 
            return self.UserK,self.UserG,self.UserRho

        
    def read_data(self):
                
        def Stinglist_float(string_list):
            number_list=[]
            
            
            def average(array):
                num=0;sum1=1e-3
                for i in array:
                    try:
                        sum1+=float(i)
                        num+=1
                    except:
                        pass
                if num == 0:
                    print (array)
                return sum1/num
                    
            
            
            def string_float(string):
                return_string = ''
                for i in string:
                    if i !=' ':
                        return_string +=i
                a=[(x) for x in re.split(r'[()]',string)]    
                list1=[]
                for i in a:
                    try:
                        list1.append(float(i))
                    except:
                        pass
                if len(list1) ==1:
                    list1.append('NA')
                if len(list1) !=2:
                    #print ('wrong input',list1)
                    list1=['NA','NA']
                return list1
            
            for i in string_list:
                list1 = string_float(i)
                number_list.append(list1)
            aa=np.array(number_list)   
            
            ave=average(aa[:,1])
            
            for i in range(len(aa[:,1])):
                if aa[:,1][i] == 'NA':
                    aa[:,1][i] = ave
            
            return aa[:,0],aa[:,1]        
        
        
        
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
           for i in range(len(line)):  
#==============================================================================
#                try:
#                    dictionary[name[i]].append(float(line[i]))
#                except:
#==============================================================================
                   dictionary[name[i]].append((line[i])) 
        file.close()
        
        
        
        self.dictionary = dictionary
        self.Flag = dictionary['Flag']
        
        a,b = Stinglist_float(dictionary['H2O'])
        self.water_content=a;
        self.water_content_error = b
        
        a,b = Stinglist_float(dictionary['Iron'])
        self.iron_content=a
        self.iron_content_error = b
        
        a,b = Stinglist_float(dictionary['K'])
        self.K = a
        self.K_error=b
        
        
        a,b = Stinglist_float(dictionary['G'])
        self.G = a
        self.G_error=b
        
        a,b = Stinglist_float(dictionary['Rho'])
        self.Rho = a
        self.Rho_error=b
        self.reference = dictionary['Author']
        self.number_of_data=len(self.K)     
        
        
        
        #print (len(self.water_content))   
        return name,dictionary  
    


    def func(self,beta,data):
        # linear function use to do the regression 
        x,y = data
        a,b,c = beta
        return a*x+b*y+c

    '''
    3 functions that return regression coefficient 
    '''
    def function_K(self,methods =2 ):
        return self.Regression_Plane(self.water_content,self.water_content_error,
                    self.iron_content,self.iron_content_error,self.K,self.K_error,methods=methods)     
    def function_G(self,methods =2):
        return self.Regression_Plane(self.water_content,self.water_content_error,
                    self.iron_content,self.iron_content_error,self.G,self.G_error,methods=methods) 
    def function_Rho(self,methods =2):
        return self.Regression_Plane(self.water_content,self.water_content_error,
                    self.iron_content,self.iron_content_error,self.Rho,self.Rho_error,methods=methods) 

    
    def Show_fit_function(self,beta,sd_beta,name,error=True):
        # function that show fit results.             
        b11=beta[0];c11=beta[1];a11=beta[2];
        b12=abs(sd_beta[0]);c12=abs(sd_beta[1]);a12=abs(sd_beta[2]);
        if error == True:
            string1=name + ' %5.1f(%5.1f)' %(a11,a12)
            def plus_Minus(b11,b12):
                if b11>=0:
                    return '+%5.1f(%5.1f)' %(b11,b12)
                else:
                    return '%5.1f(%5.1f)' %(b11,b12)
            string2=plus_Minus(b11,b12)
            string3=plus_Minus(c11,c12)
            
            string=string1+string2+'*CWater'+string3+'*XFe'
            
        else:
            string1='%.2f' %(a11)
            def plus_Minus(b11,b12):
                if b11>=0:
                    return '+%.1f' %(b11)
                else:
                    return '%.1f' %(b11)
            string2=plus_Minus(b11,b12)
            string3=plus_Minus(c11,c12)
            if name == "K'":
                string=string2+'*CWater'+string3+'*XFe+'+string1+","+name+"="+str(self.pressure_deri[0])     
            elif name == "G'":
                string=string2+'*CWater'+string3+'*XFe+'+string1+","+name+"="+str(self.pressure_deri[1])                   
            else:
                string=string2+'*CWater'+string3+'*XFe+'+string1+""+name+""+''                     
        return string 
        
    
        
    def PLOT(self,return_fig=True,methods=1):  
        '''
        This function return PLOT using data from table 
        '''
        self.regression_methods = methods
        self.Change = False
        if self.name == "Olivine":
            self.pressure_deri = [4.217,1.462]
        if self.name == 'Wadsleyite':
            self.pressure_deri = [4.322,1.444]
        if self.name == 'Ringwoodite':
            self.pressure_deri = [4.22,1.354]
        self.str1=self.Show_fit_function(self.function_K(methods = methods)[0],self.function_K(methods = methods)[1],'')
        self.str2=self.Show_fit_function(self.function_G(methods = methods)[0],self.function_G(methods = methods)[1],'')
        self.str3=self.Show_fit_function(self.function_Rho(methods = methods)[0],self.function_Rho(methods = methods)[1],'')
                
        self.d1,self.b1,self.c1,self.a1=0,self.function_K(methods = methods)[0][0],self.function_K(methods = methods)[0][1],self.function_K(methods = methods)[0][2]
        self.d2,self.b2,self.c2,self.a2=0,self.function_G(methods = methods)[0][0],self.function_G(methods = methods)[0][1],self.function_G(methods = methods)[0][2]
        self.d3,self.b3,self.c3,self.a3=0,self.function_Rho(methods = methods)[0][0],self.function_Rho(methods = methods)[0][1],self.function_Rho(methods = methods)[0][2]
        
        if return_fig ==True:        
            fig1=plot(self.water_content,self.iron_content,self.K,self.a1,self.b1,self.c1,self.d1,self.name+" Bulk modulus",self.str1)
            fig2=plot(self.water_content,self.iron_content,self.G,self.a2,self.b2,self.c2,self.d2,self.name+" Shear modulus",self.str2)
            fig3=plot(self.water_content,self.iron_content,self.Rho,self.a3,self.b3,self.c3,self.d3,self.name+" Density",self.str3)
            return fig1,fig2,fig3
        else:    
            return [self.water_content,self.iron_content,self.K,self.a1,self.b1,self.c1,self.d1,self.name+" Bulk modulus",self.str1],[self.water_content,self.iron_content,self.G,self.a2,self.b2,self.c2,self.d2,self.name+" Shear modulus",self.str2],[self.water_content,self.iron_content,self.Rho,self.a3,self.b3,self.c3,self.d3,self.name+" Density",self.str3]
     
     
    def Return_error(self):
        a = self.water_content_error
        b = self.iron_content_error
        c = self.K_error
        d = self.G_error
        e = self.Rho_error
        return [a,b,c],[a,b,d],[a,b,e]
        
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
    def __init__(self,Minerals=None,string=None,flag = None,methods = 1):
        super(Regression_PLOT_PyQt, self).__init__()
        if string is not None:
            self.stringK,self.stringG,self.stringRho,self.user_input = string
        else:
            self.stringK=None;self.stringG=None;self.stringRho=None;self.user_input=False
            
        self.regression_methods=methods
        self.resize(1400, 600)
        self.Minerals = Minerals

        self.dirty=False
        #self.user_input = user_input

        self.Minerals.original_flag()
        self.Minerals.read_data()
        self.table = QTableView()
        self.model = QStandardItemModel(self)
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
        
        self.Mehods_choice = QComboBox()
        self.Mehods_choice.addItems(['standard least squares','Orthogonal distance regression (with error)'])
        self.Mehods_choice.setCurrentIndex(self.regression_methods)
        self.Mehods_choice.currentIndexChanged.connect(self.METHODS)
        
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
            self.Kinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_K(methods =self.regression_methods)[0],self.Minerals.function_K(methods = self.regression_methods)[1],"K'",error=False))
            self.Ginput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_G(methods = self.regression_methods)[0],self.Minerals.function_G(methods = self.regression_methods)[1],"G'",error=False))
            self.Rhoinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_Rho(methods = self.regression_methods)[0],self.Minerals.function_Rho(methods = self.regression_methods)[1],'',error=False))
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
        #self.regression_methods=0
        self.PLOT()   
        self.layout.addWidget(self.table,0,1,1,17)
        self.layout.addWidget(self.button,2,0,1,9)
        self.layout.addWidget(self.button1,2,9,1,4)
        self.layout.addWidget(self.Mehods_choice,2,16,1,1)
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
  
    def METHODS(self,a):
        self.regression_methods=a
    
        
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
            self.Minerals.read_data()
            self.model.clear()
            self.model.setHorizontalHeaderLabels(['flag','Water Content','Iron Content','K (Gpa)','G (Gpa)', 'Rho (g/cm3)','Reference'])
            for i in range(len(self.Minerals.Flag)):
                a=self.Minerals.Return_original_data(i)
                item = QStandardItem(str(self.Minerals.Flag[i]))
                self.model.setItem(i, 0,item)  
                for j in range(1,7):
                    item = QStandardItem(str(a[j]))
                    self.model.setItem(i, j,item)  
                    if j != 0:
                        item.setFlags(Qt.ItemIsEnabled)
                        item.setBackground(QColor(211,211,211))
                        
        #print (self.Minerals.Change)     
            
            
        if self.user_input == False:
            self.a,self.b,self.c = self.Minerals.PLOT(return_fig=False,methods=self.regression_methods)       
            #print (self.Minerals.number_of_data)
            self.Kinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_K(methods = self.regression_methods)[0],self.Minerals.function_K(methods = self.regression_methods)[1],"K'",error=False))
            self.Ginput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_G(methods = self.regression_methods)[0],self.Minerals.function_G(methods = self.regression_methods)[1],"G'",error=False))
            self.Rhoinput_formula.setText(self.Minerals.Show_fit_function(self.Minerals.function_Rho(methods = self.regression_methods)[0],self.Minerals.function_Rho(methods = self.regression_methods)[1],'',error=False))

        else:
            self.a,self.b,self.c = self.Minerals.PLOT_input_formula(return_fig=False)   
       
        
        
        error1,error2,error3 = self.Minerals.Return_error()
            
        self.canvas1 = FigureCanvas3D(self.a,error1)
        self.canvas2 = FigureCanvas3D(self.b,error2)
        self.canvas3 = FigureCanvas3D(self.c,error3)
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
        print ('export')
        
    def ReturnString(self):
        aa = self.Minerals.Flag
        bb=str(int(aa[0]))
        for i in range(1,len(aa)):
            bb+=str(int(aa[i]))
        return [self.stringK,self.stringG,self.stringRho,self.user_input],bb,self.regression_methods



olivine=Regression('Olivine')
wadsleyte=Regression('Wadsleyite')
ringwoodite=Regression('Ringwoodite')        
        
        
if __name__ == "__main__":

    qApp = QApplication(sys.argv)
    app = Regression_PLOT_PyQt(ringwoodite,flag='')
    app.show()
    sys.exit(qApp.exec_())        
    
#==============================================================================
#     #olivine.PLOT()
#     if True:
#         Minerals = wadsleyte;cc=1;methods=2
#         Minerals.read_data()
#         a,b,c = Minerals.PLOT(return_fig=False,methods=methods)  
#         error1,error2,error3 = Minerals.Return_error()
#         if cc==1:
#             xxx = FigureCanvas3D(a,error1)
#         if cc==2:
#             xxx = FigureCanvas3D(b,error2)
#         if cc==3:
#             xxx = FigureCanvas3D(c,error3)
#         xxx.show()
# 
#==============================================================================



#print func(a[1][0],[0.1,0])#water==%,Iron<=1