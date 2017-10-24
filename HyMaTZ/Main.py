# This file is part of HyMaTZ - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2017 - 2020 by the HyMaTZ team
from __future__ import print_function
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

#import matplotlib as mpl
import  matplotlib.pyplot  as plt
from matplotlib.figure import Figure
from matplotlib import gridspec

from scipy.interpolate import UnivariateSpline
from Mineral_Physics.Extract_Data import Extrac_data
from Mineral_Physics.regression import  (Regression_PLOT_PyQt,olivine,wadsleyte,ringwoodite,Regression) 
from GUI_tools.wf_GUI_plot import (Test_Thread,wfQ_ApplicationWindow,MantleDlg)
from GUI_tools.GUI_progressbar import Progressbar
from GUI_cal_MP import Phase_diagram
from Seismic.Seismic import PREM_Suzan,AK135,IASP91,PEMC,Suzan_one,Suzan_nul,Suzan_ref
#from GUI_Phase import  Phase_extract
from GUI_Mineral_Params import Mineral_Params
from GUI_treeview import GUI_tree
from GUI_water import water_profile
from GUI_iron import iron_profile
from GUI_attenuation import Attenuation
from GUI_layer import GUI_layers,Layers

#from Mineral_Physics.Phase_diamgra1 import Olivine_phase_diagram
from Mineral_Physics.Velocity_calculator import Velocity_calculator
from Mineral_Physics.Stix2011data import (ab,an,sp,hc,fo,fa,mgwa,fewa,mgri,feri,en,fs,mgts,odi,di\
                          ,he,cen,cats,jd,hpcen,hpcfs,mgpv,fepv\
                          ,alpv,capv,mgil,feil,co,py ,al,gr,mgmj ,jdmj,qtz,coes,st,mppv,fppv,appv,pe,wu,mgcf \
                          ,fecf ,nacf,ky,neph,OL_,WA_,RI_)
from Mineral_Physics.Solidsolution import (c2c,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,ppv,Pv,Ring,Sp,Wad,OLwater,WAwater,RIwater)
from Mineral_Physics.attenuation import Anelastic
#from na04.Na04 import NA04
#from GUI_map import WorldMap

__version__ = "2.1.0"

#mpl.rcParams["figure.dpi"] = 1000
class MainWindow(QMainWindow):

    def __init__(self,parent=None):
        super(MainWindow, self).__init__(parent)        
        self.dirty = False
        self.filename = None
        self.TT = None
        self.model_choice = None
        self.model_choice_message = ''
        self.waterusercontorl={
            'OL_function':'4.0*1e-6*D*D',
            'WA_function':'3.3',
            'RI_function':'1.7',
            'OLcoefficient':0,
            'WAcoefficient':0,
            'RIcoefficient':0,
            }
        self.ironusercontorl={
            'OL_function':'88',
            'WA_function':'88',
            'RI_function':'88',

            }        
        self.dialog=None
        self.Phase_diagram = None
        self.stringOL=None
        self.stringWA=None
        self.stringRI=None
        self.dialoga1 = None
        self.dialoga2 = None
        self.dialogvalue = None
        self.olflag = '1111'
        self.waflag = '111'
        self.riflag = '111'
        self.OLOL = None
        self.WAWA = None
        self.RIRI = None
        self.olmethods=1
        self.wamethods=1
        self.rimethods=1

        self.Discontinuity410=None
        self.Discontinuity520=None
        self.Discontinuity660=None                
        #basic calculation settings
        self.sys = Phase_diagram(300) 
        self.Harzburgite = [self.sys.Depth,np.zeros(len(self.sys.Depth))]
        self.fo=np.zeros(len(self.sys.Depth))
        self.mgwa=np.zeros(len(self.sys.Depth))
        self.mgri=np.zeros(len(self.sys.Depth))
        self.ol_water=np.zeros(len(self.sys.Depth))
        self.wa_water=np.zeros(len(self.sys.Depth))
        self.ri_water=np.zeros(len(self.sys.Depth))
        self.Vp=np.zeros(len(self.sys.Depth)) 
        self.Vs=np.zeros(len(self.sys.Depth))
        #self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.plot_seismic_type_choice = 0
        self.plot_type_choice = 0
        self.average = None

        #Set title
        self.setWindowTitle("HyMaTZ")
        self.setWindowIcon(QIcon(os.path.join('GUI_tools','python.jpg')))       
        #self.main_widget = QWidget(self)
        #self.main_widget.setFocus()
        #self.setCentralWidget(self.main_widget)
        
        
        self.Layout_Plot()
        self.Layout_Control()
        self.Layout_Main()
        
        self.MenuBar()
        self.StatusBar()
        self.resize(2200,1200)
    
    '''
    Main layout 
    '''    
    def  Layout_Main(self):
        self.mainSplitter = QSplitter(Qt.Horizontal) 
        self.BTN_Splitter = QSplitter(Qt.Vertical)
        self.BTN_Splitter.addWidget(self.BTN_Widget)
        self.BTN_Splitter.addWidget(self.PLOT_Button1)
        #self.BTN_Splitter.addWidget(self.PLOT_Button2)
        
        
        self.Layout_Plot()
        self.mainSplitter.addWidget(self.BTN_Splitter)
        self.mainSplitter.addWidget(self.layout_fig)
        self.setCentralWidget(self.mainSplitter)
        
  
        self.mainSplitter.setStretchFactor(0, 1)
        self.mainSplitter.setStretchFactor(1, 3)
        self.BTN_Splitter.setStretchFactor(0, 2)
        self.BTN_Splitter.setStretchFactor(1, 2)
    '''
    Main layout end
    '''     
    
    '''
    Plot layout
    '''    
    def Layout_Plot(self):
        self.layout_fig = QSplitter(Qt.Vertical)
        self.figure = Figure(dpi=200)
        self.qmc = FigureCanvas(self.figure)        
        self.qmc.mpl_connect('button_press_event', self.OnClick)
        self.ntb = NavigationToolbar(self.qmc, self)   
        #self.ntb.save_figure
        self.layout_fig.addWidget(self.qmc)
        self.layout_fig.addWidget(self.ntb)

    '''
    Plot layout end
    '''
    
    '''
    Control layout
    '''

     
    def BTN(self):
        self.BTN_layout = QGridLayout(self)
        self.BTN_Minerals()
        self.BTN_layout.addLayout(self.BTN_Minerals_layout,0,1,4,1)
        
    def update_phase_proportion(self):
        self.sys.Phase_diagram_input()
        self.Phase_diagram = np.array(self.sys.Phase_diagram)
        self.plot_function()       

   
    def showdialog(self,message1=' ',message2=' '):
       msg = QMessageBox()
       msg.setIcon(QMessageBox.Information)    
       msg.setText(message1)
       #msg.setInformativeText("This is additional information")
       msg.setWindowTitle("MessageBox")
       msg.setDetailedText(message2)
       msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
       if msg.exec_()== QMessageBox.Yes:
           return True
       else:
           return False

           
        
    def BTN_Minerals(self):
        self.BTN_Minerals_layout = QGridLayout()  

        i=0
        i+=1  
        self.label = QLabel()
        self.label.setText("1. Select a model and click the 'Load Data for the Model' button. Alternatively, to specify a layers compositional model, check this box, click 'Load Data for the Model' and select the layers model with the pop-up window.")
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,5)
        self.Check_Layer = QCheckBox()
        self.Check_Layer.setChecked(False)
        self.Check_Layer.stateChanged.connect(self.Layer_Contotrl)
        self.Check_Layer.setLayoutDirection(Qt.RightToLeft)
        self.BTN_Minerals_layout.addWidget(self.Check_Layer,i,5,1,1)           
        
        i+=1
        self.BTN_Mantle_cambo=QComboBox()
        self.BTN_Mantle_cambo.addItem("Choose One")
        self.BTN_Mantle_cambo.addItem("Basalt(Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("Harzburgite (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("Pyrolite (Xu et al. 2008) ")
        self.BTN_Mantle_cambo.addItem("Mechanical Mixture ")
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
        self.BTN_Mantle_cambo.addItem('User input ')    
        self.BTN_Mantle_cambo.setCurrentIndex(self.BTN_Mantle_cambo.findText("Choose One"))
        self.BTN_Mantle_cambo.currentIndexChanged.connect(self.Model_select)
        self.BTN_Minerals_layout.addWidget(self.BTN_Mantle_cambo,i,0,1,6)
        
        i+=1
        self.BTN_choice = QPushButton(self)  
        self.BTN_choice.setText('Load Data for the Model')
        self.BTN_choice.clicked.connect(self.Select_Model)   
        self.BTN_Minerals_layout.addWidget(self.BTN_choice,i,0,1,6)


        i+=1
        self.label = QLabel()
        self.label.setText('2. Enter the adiabatic potential temperature at the surface, then press enter:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)
        
        i+=1
        self.textEdit = QLineEdit()
        #self.textEdit.setFixedWidth(400)
        self.textEdit.setText("")
        self.textEdit.returnPressed.connect(self.editor)
        self.BTN_Minerals_layout.addWidget(self.textEdit,i,0,1,5)       
        self.label = QLabel()
        self.label.setText('(K)')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,5,1,1)     

        
#==============================================================================
#         i+=1
#         self.Mineral_params=QPushButton()
#         self.Mineral_params.setText('View/Modify Thermoelastic Data')
#         self.Mineral_params.clicked.connect(self.mineral_params)
#         self.BTN_Minerals_layout.addWidget(self.Mineral_params,i,0,1,6)
#==============================================================================

        i+=1
        self.label = QLabel()
        self.label.resize(100,50)
        self.label.setText('3. Calculate and view the model as mineral phases versus depth :')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)

        i+=1
        self.Calculate_Button=QPushButton()
        self.Calculate_Button.setText('Show Mantle Composition versus Depth')
        self.Calculate_Button.clicked.connect(self.calculate)
        self.BTN_Minerals_layout.addWidget(self.Calculate_Button,i,0,1,6)    
        
        i+=1
        self.label = QLabel()
        self.label.setText("4. When this box is checked (the default), the elastic properties of the Mg<sub>2</sub>SiO<sub>4</sub> polymorphs will be functions of the H<sub>2</sub>O content. Uncheck this box if you want a dry mantle.")
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,5)        
        self.Check_water = QCheckBox()
        self.Check_water.setChecked(True)
        self.Check_water.stateChanged.connect(self.Use_Perplex)
        self.Check_water.setLayoutDirection(Qt.RightToLeft)
        self.BTN_Minerals_layout.addWidget(self.Check_water,i,5,1,1)        
        
        i+=1
        self.label = QLabel()
        self.label.resize(100,50)
        self.label.setText('Use this pull-down menu to enter your own function for the effect of water for a specific Mg<sub>2</sub>SiO<sub>4</sub> phases:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)
        
        i+=1
        self.BTN_NAMs_cambo=QComboBox()
        self.BTN_NAMs_cambo.addItem("Choose Mineral")
        self.BTN_NAMs_cambo.addItem("Olivine")
        self.BTN_NAMs_cambo.addItem("Wadsleyite")
        self.BTN_NAMs_cambo.addItem("Ringwoodite")
        self.BTN_NAMs_cambo.setCurrentIndex(0)
        self.BTN_NAMs_cambo.currentIndexChanged.connect(self.NAMs_select)
        self.BTN_Minerals_layout.addWidget(self.BTN_NAMs_cambo,i,0,1,6)
        
        i+=1
        self.label = QLabel()
        self.label.setText('5. Option to define your own water content profile:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)
        
        i+=1
        self.BTN_water = QPushButton()
        self.BTN_water.setText('Set Water Content Profile ')
        self.BTN_water.clicked.connect(self.Water_profile)
        self.BTN_Minerals_layout.addWidget(self.BTN_water,i,0,1,6)     


        i+=1
        self.label = QLabel()
        self.label.setText('6. Option to define your own iron content profile (for NAMs). The box must first be checked.')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,5)
        self.Check_Iron = QCheckBox()
        self.Check_Iron.setChecked(False)
        self.Check_Iron.stateChanged.connect(self.Use_Perplex)
        self.Check_Iron.setLayoutDirection(Qt.RightToLeft)        
        self.BTN_Minerals_layout.addWidget(self.Check_Iron,i,5,1,1)   
        
        i+=1
        self.BTN_iron = QPushButton()
        self.BTN_iron.setText('Set Iron Content Profile ')
        self.BTN_iron.clicked.connect(self.Iron_profile)
        self.BTN_Minerals_layout.addWidget(self.BTN_iron,i,0,1,6) 

        i+=1
        self.label = QLabel()
        self.label.resize(100,50)
        self.label.setText('7. Select averaging scheme from the next pull-down menu:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)
        

        
        i+=1
        self.AverageSchemesButton = QComboBox()
        self.AverageSchemesButton.addItem('Choose One')
        self.AverageSchemesButton.addItem('Voigt')
        self.AverageSchemesButton.addItem('Reuss')
        self.AverageSchemesButton.addItem('Voigt-Reuss-Hill')
        self.AverageSchemesButton.addItem('Hashin and Shtrikman')
        self.AverageSchemesButton.addItem('Hashin and Shtrikman upper')
        self.AverageSchemesButton.addItem('Hashin and Shtrikman lower')
        self.AverageSchemesButton.setCurrentIndex(0)
        self.AverageSchemesButton.currentIndexChanged.connect(self.average_function)
        self.BTN_Minerals_layout.addWidget(self.AverageSchemesButton,i,0,1,6)
           
        i+=1
        self.label = QLabel()
        self.label.resize(100,50)
        self.label.setText('8. Calculate Vp and Vs. Check this box to use a chosen attenuation function on the calculated velocities:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,5)   
        self.Check_attenuation = QCheckBox()
        self.Check_attenuation.setChecked(False)
        self.Check_attenuation.stateChanged.connect(self.Use_Perplex)
        self.Check_attenuation.setLayoutDirection(Qt.RightToLeft)        
        self.BTN_Minerals_layout.addWidget(self.Check_attenuation,i,5,1,1)
        

        i+=1
        self.BTN_Vp_Vs = QPushButton()
        self.BTN_Vp_Vs.setText('Calculate Vp and Vs ')
        self.BTN_Vp_Vs.clicked.connect(self.plot_calculate_velosity)
        self.BTN_Minerals_layout.addWidget(self.BTN_Vp_Vs,i,0,1,6)    

        i+=1
        self.label = QLabel()
        self.label.resize(100,50)
        self.label.setText('9. Select the plotted variable with this pull down menu:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)   
        
        i+=1
        self.Plot_type=QComboBox()
        self.Plot_type.addItem('Choose What to Plot') #0
        self.Plot_type.addItem('Clean')  #self.plot_type_choice=8
        self.Plot_type.addItem('Vp')    #self.plot_type_choice=1  
        self.Plot_type.addItem('Vs')    #self.plot_type_choice=2 
        self.Plot_type.addItem('Vp,Vs')  #self.plot_type_choice=3 
        self.Plot_type.addItem('Vp/Vs')  #self.plot_type_choice=4
        self.Plot_type.addItem('Temperature')  #self.plot_type_choice=5 
        self.Plot_type.addItem('Bulk Modulus')#self.plot_type_choice=6
        self.Plot_type.addItem('Shear Modulus')#self.plot_type_choice=7 
        self.Plot_type.addItem('Density')#self.plot_type_choice=9  
        self.Plot_type.addItem('dlnVp/dD,dlnVs/dD')
        self.Plot_type.addItem('dlnVp/dP,dlnVs/dP')
        self.Plot_type.addItem('dlnVp/dT,dlnVs/dT')
        self.Plot_type.addItem('Pressure')
        
        #self.Plot_type.addItem('dlnVp/dT,dlnVs/dWater')
        self.Plot_type.setCurrentIndex(0)
        self.Plot_typelist = ['Choose One','Vp','Vs','Vp,Vs','Vp/Vs','Temperature','Bulk Modulus','Shear Modulus','Clean','Density','dlnVp/dD,dlnVs/dD']
        self.Plot_type.currentIndexChanged.connect(self.plot_type)
        self.BTN_Minerals_layout.addWidget(self.Plot_type,i,0,1,6)
        
        i+=1
        self.label = QLabel()
        self.label.resize(100,50)
        self.label.setText('10.Add a reference model to your plot here:')
        self.label.setWordWrap(True)
        self.BTN_Minerals_layout.addWidget(self.label,i,0,1,6)
        
        i+=1
        self.Plot_seismice_type=QComboBox()
        self.Plot_seismice_type.addItem('Choose One')
        self.Plot_seismice_type.addItem('Clean')
        self.Plot_seismice_type.addItem('PREM (Dziewonski and Anderson, 1981)')
        self.Plot_seismice_type.addItem('AK135 (Kennett, Engdahl and Buland,1995) ')
        self.Plot_seismice_type.addItem('IASP91 (Kennett and Engdahl, 1991)')
        self.Plot_seismice_type.addItem('PEM-C (Dziewonski et al., 1975)' )
        self.Plot_seismice_type.addItem('User input')    
        self.Plot_seismice_type.addItem('Na04')  
        self.Plot_seismice_type.addItem('Na04-ref')      
        self.Plot_seismice_type.addItem('Na04-1') 
        self.Plot_seismice_type.addItem('Na04-2') 
        self.Plot_seismice_type.setCurrentIndex(0)
        self.Plot_seismice_type.currentIndexChanged.connect(self.plot_seismic_type)
        self.BTN_Minerals_layout.addWidget(self.Plot_seismice_type,i,0,1,6)

    def plot_type(self,a):
        if a==1:
            self.plot_type_choice=8
        if a==2:
            self.plot_type_choice=1            
        if a==3:
            self.plot_type_choice=2        
        if a==4:
            self.plot_type_choice=3  
        if a==5:
            self.plot_type_choice=4             
        if a==6:
            self.plot_type_choice=5 
        if a==7:
            self.plot_type_choice=6 
        if a==8:
            self.plot_type_choice=7           
        if a==9:
            self.plot_type_choice=9 
        if a==10:
            self.plot_type_choice=10 
        if a==11:
            self.plot_type_choice=11 
        if a==12:
            self.plot_type_choice=12 
        if a==13:
            self.plot_type_choice=13 
        if a==14:
            self.plot_type_choice=14 
            
    def average_function(self,a):
        if a==1:
            self.average=1
        if a==2:
            self.average=2        
        if a==3:
            self.average=3
        if a==4:
            self.average=4
        if a==5:
            self.average=5
        if a==6:
            self.average=6
                                
    def Layer_Contotrl(self):
        if self.Check_Layer.isChecked() == False:
            self.BTN_Mantle_cambo.setEnabled(True)
        else:
            self.BTN_Mantle_cambo.setCurrentIndex(0)
            self.BTN_Mantle_cambo.setEnabled(False)
        
    def Use_Perplex(self):
        pass
        
            
    def plot_seismic_type(self,a):
        if a==1:
            self.plot_seismic_type_choice=0
        if a==2:
            self.plot_seismic_type_choice=1
        if a==3:
            self.plot_seismic_type_choice=2 
        if a==4:
            self.plot_seismic_type_choice=3     
        if a==6:
            self.plot_seismic_type_choice=5  
            self.address_input_seismic, _  = QFileDialog.getOpenFileName(self, 'Open file', '/home')
        if a==5:
            self.plot_seismic_type_choice=4
        if a ==7:
            self.plot_seismic_type = 6
            try:
                self.dialog = WorldMap(self.na04)
            except:
                self.na04  = NA04()
                self.dialog = WorldMap(self.na04)
            if self.dialog.exec_():
                self.address_input_seismic = os.path.join(os.getcwd(),'result.txt')                
                    
        if a==8:
            self.plot_seismic_type_choice=7                
        if a==9:
            self.plot_seismic_type_choice=8          
        if a==10:
            self.plot_seismic_type_choice=9

    def Iron_profile(self):
        if self.Check_Iron.isChecked() != True:
            pass
        else:
            self.dialog1= iron_profile(Depth =self.sys.Depth,Pressure =self.sys.Pressure, Temperature = self.sys.Temperature,usercontrol=self.ironusercontorl,fo=self.fo,mgwa=self.mgwa,mgri=self.mgri)
            if self.dialog1.exec_():
                self.ironusercontorl = self.dialog1.usercontrol    
            self.fo = self.dialog1.OL_iron
            self.mgwa = self.dialog1.WA_iron
            self.mgri = self.dialog1.RI_iron

        
    def Water_profile(self):
        if self.Check_water.isChecked() != True:
            pass
        else:
            self.dialog= water_profile(Depth =self.sys.Depth,Pressure =self.sys.Pressure, Temperature = self.sys.Temperature,usercontrol=self.waterusercontorl)
            if self.dialog.exec_():
                self.waterusercontorl = self.dialog.usercontrol

    def mineral_params(self):
        dialog= Mineral_Params(VolumeOff=True)
        if dialog.exec_():
            pass

    
    def NAMs_select(self,a):
 
        if self.Check_water.isChecked() != True:
            pass
        else:   
            if a == 1:
                self.appOL = Regression_PLOT_PyQt(olivine,self.stringOL,self.olflag,self.olmethods)
                if self.appOL.exec_():
                    pass
                self.stringOL,self.olflag,self.olmethods = self.appOL.ReturnString()
                #print (self.olflag)

            if a == 2:
                self.appWA = Regression_PLOT_PyQt(wadsleyte,self.stringWA,self.waflag,self.wamethods)
                if self.appWA.exec_():
                    pass
                self.stringWA,self.waflag,self.wamethods= self.appWA.ReturnString()   
                #print (self.waflag)

            if a == 3:
                self.appRI = Regression_PLOT_PyQt(ringwoodite,self.stringRI,self.riflag,self.rimethods)
                if self.appRI.exec_():
                    pass
                self.stringRI,self.riflag,self.rimethods = self.appRI.ReturnString()
            self.BTN_NAMs_cambo.setCurrentIndex(0)  
            #self.stringOL = self.appOL.ReturnString()         
            
            #self.stringOL=[self.appOL.stringK,self.appOL.stringG,self.appOL.stringRho]
            #self.stringWA=[self.appWA.stringK,self.appWA.stringG,self.appWA.stringRho]
            #self.stringRI=[self.appRI.stringK,self.appRI.stringG,self.appRI.stringRho]


        

        
    def Model_select(self,a):  
        #self.message=''
        list1 = ['CaO', 'Al2O3',  'Na2O', 'MgO', 'FeO', 'SiO2']
        list2 = np.array([  56.0774,  101.9612,   61.979 ,   40.3044,   71.8444,   60.0843])
        Harzburgite = np.array([0.81, 0.53, 0.00, 56.51, 6.07, 36.07])
        Basalt      = np.array([13.88, 10.19, 2.18, 14.94, 7.06, 51.75])
        Pyrolite    = np.array([2.93553,2.21491,0.109863,49.8509,6.16912,38.6950])
        if a == 1:
            self.Type = 'Basalt'
            self.message = 'Basalt (Xu et al. 2008) '
            self.savemessage = 'Basalt (Xu et al. 2008) '
            for i in range(len(self.Harzburgite[0])):
                self.Harzburgite[1][i] = 0
            self.Model_composition = Basalt 
        if a == 2:
            self.Type = 'Harzburgite100'
            self.message = 'Harzburgite (Xu et al. 2008) '
            self.savemessage = 'Harzburgite (Xu et al. 2008) '
            for i in range(len(self.Harzburgite[0])):
                self.Harzburgite[1][i] = 1
            self.Model_composition = Harzburgite
        if a == 3:
            self.Type = 'Pyrolite'
            self.message = 'Pyrolite (Xu et al. 2008) '
            self.savemessage = 'Pyrolite (Xu et al. 2008) '
            for i in range(len(self.Harzburgite[0])):
                self.Harzburgite[1][i] = 1.01  
            self.Model_composition =  Pyrolite      
        if a == 4:
            self.Type = 'Mechanical Mixture'
            self.message = 'Mechanical Mixture'
            self.savemessage = 'Mechanical Mixture'
            
            self.dialog = MantleDlg(self.dialoga1,self.dialoga2,self.dialogvalue)
            if self.dialog.exec_():
                pass
            self.dialoga1 = self.dialog.a1
            self.dialoga2 = self.dialog.a2
            self.dialogvalue = self.dialog.value
            self.Model_composition =  copy.deepcopy(self.dialog.composition)
            #print ("Done")
                
                #print (aa)
        if a == 5:
            self.Type = 'Piclogite'
            self.message = 'Piclogite (Weidner, 1986) '
            self.savemessage = 'Piclogite (Weidner, 1986) '
            self.Model_composition = np.array([  7.22,8.16,2.80,30.14,4.58,47.07])
            for i in range(len(self.Harzburgite[0])):
                self.Harzburgite[1][i] = 1.02   
        if a ==6:
            self.Type = 'Harzburgite90'
            self.message = '90%molHarzburgite + 10%molBasalt (Xu et al. 2008) '
            self.savemessage = '90%molHarzburgite + 10%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.9
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt            
        if a ==7:
            self.Type = 'Harzburgite80'
            self.message = '80%molHarzburgite + 20%molBasalt (Xu et al. 2008) '
            self.savemessage = '80%molHarzburgite + 20%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.8
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==8:
            self.Type = 'Harzburgite70'
            self.message = '70%molHarzburgite + 30%molBasalt (Xu et al. 2008) '
            self.savemessage = '70%molHarzburgite + 30%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.7
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==9:
            self.Type = 'Harzburgite60'
            self.message = '60%molHarzburgite + 40%molBasalt (Xu et al. 2008) '
            self.savemessage = '60%molHarzburgite + 40%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.6
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==10:
            self.Type = 'Harzburgite50'
            self.message = '50%molHarzburgite + 50%molBasalt (Xu et al. 2008) '
            self.savemessage = '50%molHarzburgite + 50%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.5
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==11:
            self.Type = 'Harzburgite40'
            self.message = '40%molHarzburgite + 60%molBasalt (Xu et al. 2008) '
            self.savemessage = '40%molHarzburgite + 60%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.4
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==12:
            self.Type = 'Harzburgite30'
            self.message = '30%molHarzburgite + 70%molBasalt (Xu et al. 2008) '
            self.savemessage = '30%molHarzburgite + 70%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.3
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==13:
            self.Type = 'Harzburgite20'
            self.message = '20%molHarzburgite + 80%molBasalt (Xu et al. 2008) '
            self.savemessage = '20%molHarzburgite + 80%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.2
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt  
        if a ==14:
            self.Type = 'Harzburgite10'
            self.message = '10%molHarzburgite + 90%molBasalt (Xu et al. 2008) '
            self.savemessage = '10%molHarzburgite + 90%molBasalt (Xu et al. 2008) '
            Harzburgite_precentage = 0.1
            self.Model_composition = Harzburgite_precentage*Harzburgite + (1-Harzburgite_precentage)*Basalt 
        if a ==15:
            self.Type = 'Pyrolite1'
            self.message = 'Pyrolite1 (McDonough and Sun, 1995) '
            self.savemessage = 'Pyrolite1 (McDonough and Sun, 1995) '
            self.Model_composition = np.array([  3.31,2.28,0.30,49.06,5.86,39.18])
        if a ==16:
            self.Type = 'Pyrolite2'
            self.message = 'Pyrolite2 (Jagoutz et al., 1979) '
            self.savemessage = 'Pyrolite2 (Jagoutz et al., 1979) '
            self.Model_composition = np.array([  3.25,2.03,0.27,49.56,5.67,39.19])
        if a ==17:
            self.Type = 'Pyrolite3'
            self.message = 'Pyrolite3 (Green et al., 1979) '
            self.savemessage = 'Pyrolite3 (Green et al., 1979) '
            self.Model_composition = np.array([  3.14,2.23,0.33,49.94,5.48,38.85])
        if a ==18:
            self.Type = 'Pyrolite4'
            self.message = 'Pyrolite4 (Taylor and McLennan, 1985) '
            self.savemessage = 'Pyrolite4 (Taylor and McLennan, 1985) '
            self.Model_composition = np.array([  2.71,1.87,0.28, 45.73,5.83,43.55])

        if a ==19:
            self.Type = 'Pyrolite5'
            self.message = 'Pyrolite5 (Ringwoode, 1962) '
            self.savemessage = 'Pyrolite5 (Ringwoode, 1962) '
            self.Model_composition = np.array([  2.90,1.82,0.50, 48.96,6.22,39.58])
            
        if a ==20:
            self.Type = 'Input'
            self.message = 'User Input '
            self.savemessage = 'User Input '
            self.Model_composition = np.array([  2.90,1.82,0.50, 48.96,6.22,39.58])  
        
#==============================================================================
#         if a ==21:
#             self.Type = 'Mars'
#             self.message = 'Mars'
#             self.savemessage = 'Mars'
#             self.Model_composition = np.array([ 0,0,0,0,0,0]) 
# 
#         if a ==22:
#             self.Type = 'PyrMg-2'
#             self.message = 'PyrMg-2'
#             self.savemessage = 'PyrMg-2'
#             self.Model_composition = np.array([ 0,0,0,0,0,0]) 
# 
#         if a ==23:
#             self.Type = 'PyrMg-1'
#             self.message = 'PyrMg-1'
#             self.savemessage = 'PyrMg-1'
#             self.Model_composition = np.array([ 0,0,0,0,0,0]) 
#             
#         if a ==24:
#             self.Type = 'Pyr'
#             self.message = 'Pyr'
#             self.savemessage = 'Pyr'
#             self.Model_composition = np.array([ 0,0,0,0,0,0])    
#             
#         if a ==25:
#             self.Type = 'PyrMg+1'
#             self.message = 'PyrMg+1'
#             self.savemessage = 'PyrMg+1'
#             self.Model_composition = np.array([ 0,0,0,0,0,0])
# 
#         if a ==26:
#             self.Type = 'PyrMg+2'
#             self.message = 'PyrMg+2'
#             self.savemessage = 'PyrMg+2'
#             self.Model_composition = np.array([ 0,0,0,0,0,0])
# 
#         if a ==27:
#             self.Type = 'PyrMg+3'
#             self.message = 'PyrMg+3'
#             self.savemessage = 'PyrMg+3'
#             self.Model_composition = np.array([ 0,0,0,0,0,0])
#==============================================================================

             
        self.Model_weight_composition = self.Model_composition * list2
        sum1 = sum(self.Model_weight_composition)
        self.Model_weight_composition /= sum1/100
        
        for i in range(6):           
            self.message += list1[i]
            self.message += ' ['
            self.message += "{0:.2f}".format((self.Model_composition[i]))
            self.message += '/'
            self.message += "{0:.2f}".format((self.Model_weight_composition[i]))
            self.message += ']; '
        self.message += 'Bulk composition in mol%/wt%'
        self.status.showMessage(self.message)  
        self.model_choice=self.Type
        return self.Type             

        
    def Select_Model(self):     
        if self.Check_Layer.isChecked() is False:
            self.BTN_Mantle_cambo.setEnabled(True)
            self.BTN_Mantle_cambo.setCurrentIndex(0)
            self.progress = QProgressDialog("Running","Cancel",0,30,self) 
            self.progress.setWindowTitle('Please wait...')
            self.progress.setWindowModality(Qt.WindowModal)
            self.progress.setMinimumDuration(0)
            self.progress.setValue(0)
            self.progress.canceled.connect(self.progress.close)
            for i in range(15):
                self.progress.setValue(i)
            self.progress.deleteLater()
            if self.Type == 'Input':
                address_input, _  = QFileDialog.getOpenFileName(self, 'Open file', '/home')
                address_output2 = address_input[:-5]+'_2.txt'
                address_output3 = address_input[:-5]+'_3.txt'                
                Extrac_data(address_input,address_output2,address_output3)  
            
            if self.model_choice == 'Mechanical Mixture':
                try:
                    aa = self.dialog.value/100.
                    a1 = self.dialog.Type
                    a2 = self.dialog.Type1
                    #print (a1,a2,aa)
                    self.Type = [a1,a2,aa]
                    self.model_choice='Mechanical Mixture'
                except:
                    a1,a2,aa = self.Type
                self.sys.change_model1(a1,a2,aa,index_P=600,index_T=400)
                #print ('done')
            else:
                self.sys.change_model(Model=self.Type,Model_composition=self.Model_composition,index_P=600,index_T=400)
            self.model_choice_message=copy.deepcopy(self.message)
            for i in range(16,30):
                self.progress.setValue(i)
            self.progress.close()
        else:
            try:
                layer_depth,layer_name,layer_number = self.LayersDialog.Return()
                self.LayersDialog = GUI_layers(layer_depth,layer_name,layer_number)
            except:
                self.LayersDialog = GUI_layers()
            if self.LayersDialog.exec_():
                pass
            self.BTN_Mantle_cambo.setEnabled(False) 
            layer_depth,layer_name,layer_number = self.LayersDialog.Return()
            self.progress = QProgressDialog("Running","Cancel",0,30,self) 
            self.progress.setWindowTitle('Please wait...')
            self.progress.setWindowModality(Qt.WindowModal)
            self.progress.setMinimumDuration(0)
            self.progress.setValue(0)
            self.progress.canceled.connect(self.progress.close)
            for i in range(15):
                self.progress.setValue(i)
            self.progress.deleteLater()
            if self.LayersDialog.Change == True:           
                self.layers = Layers(layer_name,layer_depth)  
            else:
                pass
            for i in range(16,30):
                self.progress.setValue(i)
            self.progress.close()    
            self.model_choice ='Layer'
            self.sys.Import_layer(self.layers)
            





     
    def calculate(self,Pass=False):
        
        if self.TT is None:
            self.showdialog(message1="Enter the temperature, then press enter",message2="Enter the adiabatic potential temperature at the surface, then press enter and click Show/Update Phase Proportion")
            return None
        if self.model_choice is None:
            self.showdialog(message1="Please select a model",message2="After select model please press select model") 
            return None        
        if Pass == False:
            if self.Check_Layer.isChecked() == True:
                self.model_choice = 'Layers model'; self.model_choice_message = 'you are using Layers model'
            if self.showdialog(message1="You are using " + self.model_choice +' model',message2=self.model_choice_message+' foot temperature:'+str(self.TT)+'(K)'):
                #cc = self.Harzburgite
                        
                self.progress = QProgressDialog("Running","Cancel",0,100,self) 
                self.progress.setWindowTitle('Please wait...')
                self.progress.setWindowModality(Qt.WindowModal)
                self.progress.setMinimumDuration(0)
                self.progress.setValue(0)
                self.progress.canceled.connect(self.progress.close)
                for i in range(10):
                    self.progress.setValue(i)
                
                #self.sys.composition_profile(composition=cc)
                if self.Check_Layer.isChecked() == True:
                    self.sys.Phase_diagram_first_pinciple(self.TT,layers=True)   
                else:
                    self.sys.Phase_diagram_first_pinciple(self.TT)       
                #print ('calculate')
                
                for i in range(10,50):
                    self.progress.setValue(i)
                    
                for number,phase in enumerate(self.sys.Phase_diagram):
                    self.fo[number]=phase[24]/(phase[24]+phase[25]+1e-11);
                    self.mgwa[number]=phase[26]/(phase[26]+phase[27]+1e-11);
                    self.mgri[number]=phase[28]/(phase[28]+phase[29]+1e-11);
                #print ('wf')
                
                for i in range(50,90):
                    self.progress.setValue(i)
                    
                self.Phase_diagram = np.array(self.sys.Phase_diagram)
                self.plot_function()
                for i in range(90,101):
                    self.progress.setValue(i)            
                #self.dialog= water_profile(Depth =self.sys.Depth,Pressure =self.sys.Pressure, Temperature = self.sys.Temperature)
        else:
                self.progress = QProgressDialog("Running","Cancel",0,100,self) 
                self.progress.setWindowTitle('Please wait...')
                self.progress.setWindowModality(Qt.WindowModal)
                self.progress.setMinimumDuration(0)
                self.progress.setValue(0)
                self.progress.canceled.connect(self.progress.close)
                for i in range(10):
                    self.progress.setValue(i)
                
                #self.sys.composition_profile(composition=cc)
                if self.Check_Layer.isChecked() == True:
                    self.sys.Phase_diagram_first_pinciple(self.TT,layers=True)   
                else:
                    self.sys.Phase_diagram_first_pinciple(self.TT)          
                #print ('calculate')
                
                for i in range(10,50):
                    self.progress.setValue(i)
                    
                for number,phase in enumerate(self.sys.Phase_diagram):
                    self.fo[number]=phase[24]/(phase[24]+phase[25]+1e-11);
                    self.mgwa[number]=phase[26]/(phase[26]+phase[27]+1e-11);
                    self.mgri[number]=phase[28]/(phase[28]+phase[29]+1e-11);
                #print ('wf')
                
                for i in range(50,90):
                    self.progress.setValue(i)
                    
                self.Phase_diagram = np.array(self.sys.Phase_diagram)
                self.plot_function()
                for i in range(90,101):
                    self.progress.setValue(i)                
            
            
            
            
         
    def plot_calculate_velosity(self):
        if self.Phase_diagram is None:
            self.showdialog(message1="Click shows phase proportion",message2="Click 'shows phase proportion' and check if enter temperature and select a model")
            return None            
        nnumber = len(self.Phase_diagram)
        if self.average is None:
            self.showdialog(message1="Choose average schemes",message2="please choose an average schemes")
            return None
        if self.dialog is None and self.Check_water.isChecked() == True:
            self.showdialog(message1="Please set up water profile",message2="Please set up water profile")
            return None            
        
        if self.Check_attenuation.isChecked() == True:
            try:
                self.attenua = Attenuation(parameter_list=self.AttenuationParameters)
            except:
                self.attenua = Attenuation()
            if self.attenua.exec_():
                pass
        
        
        dp=1e5;dt=5
        self.progress = QProgressDialog("Running","Cancel",0,nnumber-1,self) 
        self.progress.setWindowTitle('Please wait...')
        self.progress.setWindowModality(Qt.WindowModal)
        self.progress.setMinimumDuration(0)
        self.progress.setValue(0)
        self.progress.canceled.connect(self.progress.close)
        #self.progress.deleteLater()
        
        self.Vp = np.zeros(len(self.Vp));self.Vs = np.zeros(len(self.Vp));
        self.K=np.zeros(len(self.Vp));self.G=np.zeros(len(self.Vp));self.Rho=np.zeros(len(self.Vp))
        
        self.VpPU = np.zeros(len(self.Vp));self.VsPU = np.zeros(len(self.Vp));
        self.KPU=np.zeros(len(self.Vp));self.GPU=np.zeros(len(self.Vp));self.RhoPU=np.zeros(len(self.Vp))        
        self.VpPD = np.zeros(len(self.Vp));self.VsPD = np.zeros(len(self.Vp));
        self.KPD=np.zeros(len(self.Vp));self.GPD=np.zeros(len(self.Vp));self.RhoPD=np.zeros(len(self.Vp))

        self.VpTU = np.zeros(len(self.Vp));self.VsTU = np.zeros(len(self.Vp));
        self.KTU=np.zeros(len(self.Vp));self.GTU=np.zeros(len(self.Vp));self.RhoTU=np.zeros(len(self.Vp))        
        self.VpTD = np.zeros(len(self.Vp));self.VsTD = np.zeros(len(self.Vp));
        self.KTD=np.zeros(len(self.Vp));self.GTD=np.zeros(len(self.Vp));self.RhoTD=np.zeros(len(self.Vp))
        
        if self.Check_water.isChecked() != True:
            mineral_list=[an,ab,sp,hc,en,fs,mgts,odi,hpcen,hpcfs,di,he,cen,cats,jd,py,al,gr,mgmj,jdmj, \
                          capv,fo,fa,mgwa,fewa,mgri,feri,mgil,feil,co,mgpv,fepv,alpv,mppv,fppv,appv,mgcf,fecf, \
                          nacf,pe,wu,qtz,coes,st,an,ky,neph]     
            for i in mineral_list:
                i.Clear_Vp_Vs()
            for number,phase in enumerate(self.Phase_diagram):

                #phase[53]=0;phase[54]=0;
                K=[];G=[];Rho=[];V=[] 
                KPU=[];GPU=[];RhoPU=[]
                KPD=[];GPD=[];RhoPD=[]
                KTU=[];GTU=[];RhoTU=[] 
                KTD=[];GTD=[];RhoTD=[]
                
                for num,i in enumerate(mineral_list,start=3):
                    if phase[num] != 0:    
                        r,k,g,v,rho = i.EOS(self.sys.Pressure[number]*1e5,self.sys.Temperature[number])
                        K.append(k);G.append(g);Rho.append(rho);V.append(phase[num])
                        
                        rPU,kPU,gPU,vPU,rhoPU = i.EOS(self.sys.Pressure[number]*1e5+dp,self.sys.Temperature[number])
                        KPU.append(kPU);GPU.append(gPU);RhoPU.append(rhoPU)                      
                        rPD,kPD,gPD,vPD,rhoPD = i.EOS(self.sys.Pressure[number]*1e5-dp,self.sys.Temperature[number])
                        KPD.append(kPD);GPD.append(gPD);RhoPD.append(rhoPD)                             
                        
                        rTU,kTU,gTU,vTU,rhoTU = i.EOS(self.sys.Pressure[number]*1e5,self.sys.Temperature[number]+dt)
                        KTU.append(kTU);GTU.append(gTU);RhoTU.append(rhoTU)                      
                        rTD,kTD,gTD,vTD,rhoTD = i.EOS(self.sys.Pressure[number]*1e5,self.sys.Temperature[number]-dt)
                        KTD.append(kTD);GTD.append(gTD);RhoTD.append(rhoTD)   

                        a,b,c=i.Vp_Vs(self.sys.Pressure[number]*1e5,self.sys.Temperature[number])   
                        c*=1000
                        i.Store_Vp_Vs(a,b,c,phase[1],phase[num])
                            
                self.progress.setValue(number)
                KGRho = Velocity_calculator(K,G,Rho,V)
                if self.average == 1:
                    self.K[number],self.G[number],self.Rho[number] = KGRho.Voigt()
                if self.average == 2:
                    self.K[number],self.G[number],self.Rho[number] = KGRho.Reuss()
                if self.average == 3:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()  
                self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.
                #print ('wf')
                KGRhoPU = Velocity_calculator(KPU,GPU,RhoPU,V)
                if self.average == 1:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] = KGRhoPU.Voigt()
                if self.average == 2:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] = KGRhoPU.Reuss()
                if self.average == 3:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Hashinand_Shtrikman_lower()  
                self.VpPU[number],self.VsPU[number]=np.sqrt((self.KPU[number]+4.*self.GPU[number]/3.)/self.RhoPU[number])/1000. , np.sqrt(self.GPU[number]/self.RhoPU[number])/1000.
                #print ('wf')
#==============================================================================
#                 KGRhoPD = Velocity_calculator(KPD,GPD,RhoPD,V)
#                 if self.average == 1:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] = KGRhoPD.Voigt()
#                 if self.average == 2:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] = KGRhoPD.Reuss()
#                 if self.average == 3:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] =  KGRhoPD.Voigt_Reuss_Hill() 
#                 if self.average == 4:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] =  KGRhoPD.Hashinand_Shtrikman() 
#                 if self.average == 5:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] =  KGRhoPD.Hashinand_Shtrikman_upper() 
#                 if self.average == 6:
#                     self.KPD[number],self.GPD[number],self.RhoPU[number] =  KGRhoPD.Hashinand_Shtrikman_lower()  
#                 self.VpPD[number],self.VsPD[number]=np.sqrt((self.KPD[number]+4.*self.GPD[number]/3.)/self.RhoPD[number])/1000. , np.sqrt(self.GPD[number]/self.RhoPD[number])/1000.
#==============================================================================
                KGRhoTU = Velocity_calculator(KTU,GTU,RhoTU,V)
                if self.average == 1:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] = KGRhoTU.Voigt()
                if self.average == 2:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] = KGRhoTU.Reuss()
                if self.average == 3:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Hashinand_Shtrikman_lower()  
                self.VpTU[number],self.VsTU[number]=np.sqrt((self.KTU[number]+4.*self.GTU[number]/3.)/self.RhoTU[number])/1000. , np.sqrt(self.GTU[number]/self.RhoTU[number])/1000.
                #print ('wf')
                #print (KTD,GTD,RhoTD)
                KGRhoTD = Velocity_calculator(KTD,GTD,RhoTD,V)
                if self.average == 1:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] = KGRhoTD.Voigt()
                if self.average == 2:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] = KGRhoTD.Reuss()
                if self.average == 3:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Hashinand_Shtrikman_lower()  
                self.VpTD[number],self.VsTD[number]=np.sqrt((self.KTD[number]+4.*self.GTD[number]/3.)/self.RhoTD[number])/1000. , np.sqrt(self.GTD[number]/self.RhoTD[number])/1000.
                #print ('wf')

            #self.dlnVsWater = np.zeros(self.Vp);self.dlnVsWater = np.zeros(self.Vp)
        else:
            OL=OL_(olivine)
            WA=WA_(wadsleyte)
            RI=RI_(ringwoodite) 
            if self.Check_Iron.isChecked() != True:
                for number,phase in enumerate(self.sys.Phase_diagram):
                    self.fo[number]=phase[24]/(phase[24]+phase[25]+1e-11);
                    self.mgwa[number]=phase[26]/(phase[26]+phase[27]+1e-11);
                    self.mgri[number]=phase[28]/(phase[28]+phase[29]+1e-11);
            mineral_list=[an,ab,sp,hc,en,fs,mgts,odi,hpcen,hpcfs,di,he,cen,cats,jd,py,al,gr,mgmj,jdmj, \
                          capv,fo,fa,mgwa,fewa,mgri,feri,mgil,feil,co,mgpv,fepv,alpv,mppv,fppv,appv,mgcf,fecf, \
                          nacf,pe,wu,qtz,coes,st,an,ky,neph,OL,WA,RI]        
            for i in mineral_list:
                i.Clear_Vp_Vs()
            
            try:
                a = self.dialog.OL_water[0]
            except:
                self.showdialog(message1="Please set up water profile",message2="Please set up water profile")
                return None    

            for number,phase in enumerate(self.Phase_diagram):

                #OL.Set_Water_Iron_Condition(self.ol_water[number],self.fo[num])
                phase[53]=0;phase[54]=0;
                #phase[24]=0;phase[25]=0
                #phase[26]=0;phase[27]=0
               # phase[28]=0;phase[29]=0

                    
                OL.Set_Water_Iron_Condition(self.dialog.OL_water[number],1-self.fo[number])
                WA.Set_Water_Iron_Condition(self.dialog.WA_water[number],1-self.mgwa[number])
                RI.Set_Water_Iron_Condition(self.dialog.RI_water[number],1-self.mgri[number]) 
                
                K=[];G=[];Rho=[];V=[] 
                KPU=[];GPU=[];RhoPU=[]
                KPD=[];GPD=[];RhoPD=[]
                KTU=[];GTU=[];RhoTU=[] 
                KTD=[];GTD=[];RhoTD=[]

                
                for num,i in enumerate(mineral_list,start=3):
                    if phase[num] != 0:   
                        if num==24 or num==25 or num==26 or num==27 or num==28 or num==29:
                            pass
                        else:
                            r,k,g,v,rho = i.EOS(self.sys.Pressure[number]*1e5,self.sys.Temperature[number])
                            K.append(k);G.append(g);Rho.append(rho);V.append(phase[num])
                            
                            rPU,kPU,gPU,vPU,rhoPU = i.EOS(self.sys.Pressure[number]*1e5+dp,self.sys.Temperature[number])
                            KPU.append(kPU);GPU.append(gPU);RhoPU.append(rhoPU)                      
                            rPD,kPD,gPD,vPD,rhoPD = i.EOS(self.sys.Pressure[number]*1e5-dp,self.sys.Temperature[number])
                            KPD.append(kPD);GPD.append(gPD);RhoPD.append(rhoPD)                             
                            
                            rTU,kTU,gTU,vTU,rhoTU = i.EOS(self.sys.Pressure[number]*1e5,self.sys.Temperature[number]+dt)
                            KTU.append(kTU);GTU.append(gTU);RhoTU.append(rhoTU)                      
                            rTD,kTD,gTD,vTD,rhoTD = i.EOS(self.sys.Pressure[number]*1e5,self.sys.Temperature[number]-dt)
                            KTD.append(kTD);GTD.append(gTD);RhoTD.append(rhoTD)  
     
                            
                            a,b,c=i.Vp_Vs(self.sys.Pressure[number]*1e5,self.sys.Temperature[number]) 
                            c*=1000
                            i.Store_Vp_Vs(a,b,c,phase[1],phase[num])
                self.progress.setValue(number)
                
                #self.Vp[number],self.Vs[number]=phase[53],phase[54] 
                self.progress.setValue(number)
                KGRho = Velocity_calculator(K,G,Rho,V)
                if self.average == 1:
                    self.K[number],self.G[number],self.Rho[number] = KGRho.Voigt()
                if self.average == 2:
                    self.K[number],self.G[number],self.Rho[number] = KGRho.Reuss()
                if self.average == 3:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()  
                self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.
                #print ('wf')
                KGRhoPU = Velocity_calculator(KPU,GPU,RhoPU,V)
                if self.average == 1:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] = KGRhoPU.Voigt()
                if self.average == 2:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] = KGRhoPU.Reuss()
                if self.average == 3:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.KPU[number],self.GPU[number],self.RhoPU[number] =  KGRhoPU.Hashinand_Shtrikman_lower()  
                self.VpPU[number],self.VsPU[number]=np.sqrt((self.KPU[number]+4.*self.GPU[number]/3.)/self.RhoPU[number])/1000. , np.sqrt(self.GPU[number]/self.RhoPU[number])/1000.
                #print ('wf')
#==============================================================================
#                 KGRhoPD = Velocity_calculator(KPD,GPD,RhoPD,V)
#                 if self.average == 1:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] = KGRhoPD.Voigt()
#                 if self.average == 2:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] = KGRhoPD.Reuss()
#                 if self.average == 3:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] =  KGRhoPD.Voigt_Reuss_Hill() 
#                 if self.average == 4:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] =  KGRhoPD.Hashinand_Shtrikman() 
#                 if self.average == 5:
#                     self.KPD[number],self.GPD[number],self.RhoPD[number] =  KGRhoPD.Hashinand_Shtrikman_upper() 
#                 if self.average == 6:
#                     self.KPD[number],self.GPD[number],self.RhoPU[number] =  KGRhoPD.Hashinand_Shtrikman_lower()  
#                 self.VpPD[number],self.VsPD[number]=np.sqrt((self.KPD[number]+4.*self.GPD[number]/3.)/self.RhoPD[number])/1000. , np.sqrt(self.GPD[number]/self.RhoPD[number])/1000.
#==============================================================================
                KGRhoTU = Velocity_calculator(KTU,GTU,RhoTU,V)
                if self.average == 1:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] = KGRhoTU.Voigt()
                if self.average == 2:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] = KGRhoTU.Reuss()
                if self.average == 3:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.KTU[number],self.GTU[number],self.RhoTU[number] =  KGRhoTU.Hashinand_Shtrikman_lower()  
                self.VpTU[number],self.VsTU[number]=np.sqrt((self.KTU[number]+4.*self.GTU[number]/3.)/self.RhoTU[number])/1000. , np.sqrt(self.GTU[number]/self.RhoTU[number])/1000.
                #print ('wf')
                #print (KTD,GTD,RhoTD)
                KGRhoTD = Velocity_calculator(KTD,GTD,RhoTD,V)
                if self.average == 1:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] = KGRhoTD.Voigt()
                if self.average == 2:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] = KGRhoTD.Reuss()
                if self.average == 3:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Voigt_Reuss_Hill() 
                if self.average == 4:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Hashinand_Shtrikman() 
                if self.average == 5:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Hashinand_Shtrikman_upper() 
                if self.average == 6:
                    self.KTD[number],self.GTD[number],self.RhoTD[number] =  KGRhoTD.Hashinand_Shtrikman_lower()  
                self.VpTD[number],self.VsTD[number]=np.sqrt((self.KTD[number]+4.*self.GTD[number]/3.)/self.RhoTD[number])/1000. , np.sqrt(self.GTD[number]/self.RhoTD[number])/1000.
                #print ('wf')
            self.OLOL = OL
            self.WAWA = WA
            self.RIRI = RI        
        
        self.lnVp = np.log(self.Vp);self.lnVpP = np.log(self.VpPU)-np.log(self.lnVp);self.lnVpT = np.log(self.VpTU)-np.log(self.VpTD);#self.lnVpPD = np.log(self.VpPD)
        self.lnVs = np.log(self.Vs);self.lnVsP = np.log(self.VsPU)-np.log(self.lnVs);self.lnVsT = np.log(self.VsTU)- np.log(self.VsTD);#self.lnVsPD = np.log(self.VsPD)
        self.dD = self.sys.Depth[10]-self.sys.Depth[9]
        self.dlnVp = np.gradient(self.lnVp)/self.dD
        self.dlnVs = np.gradient(self.lnVs)/self.dD
        self.dlnVpP = np.gradient(self.lnVpP)/(2*dp);self.dlnVpT = np.gradient(self.lnVpT)/(2*dt);
        self.dlnVsP = np.gradient(self.lnVsP)/(2*dp);self.dlnVsT = np.gradient(self.lnVsT)/(2*dt);
        
        if self.Check_attenuation.isChecked() == True:
            if self.attenua.Type == 1:
                #print ('attenuation started')
                A=self.attenua.A
                w=self.attenua.w
                a=self.attenua.a
                H=self.attenua.H
                V=self.attenua.V
                self.AttenuationParameters = [A,w,a,H,V,0,0]
                Qs,Qp,a,b = Anelastic.Anelastic1(P=self.sys.Pressure,T=self.sys.Temperature,A=A,w=w,a=a,H=H,V=V,Vp=self.Vp,Vs=self.Vs) 
                self.Vs *= a
                self.Vp *= b
                #print ('attenuation done')
            if self.attenua.Type == 0:
                #print ('attenuation started')
                A=self.attenua.A
                w=self.attenua.w
                a=self.attenua.a
                g1=self.attenua.g1
                g2=self.attenua.g2
                self.AttenuationParameters = [A,w,a,0,0,g1,g2]
                Qs,Qp,a,b = Anelastic.Anelastic(P=self.sys.Pressure,T=self.sys.Temperature,B=A,w=w,a=a,g1=g1,g2=g2,Vp=self.Vp,Vs=self.Vs) 
                self.Vs *= a
                self.Vp *= b
                #print ('attenuation done')

        
        
     
    
    def editor(self):
        self.TT = (float(self.textEdit.text())) 
        if self.TT <=800 or self.TT >=2000:
            self.showdialog('re-enter temperature','temperature should between 800-2000K, here temperature will change to 1600k')
            self.textEdit.setText('1600')
            self.TT = 1600
        #self.sys.Phase_diagram_first_pinciple(self.TT)

        
    def setMantleProperties(self):
        dialog = MantleDlg(self)
        if dialog.exec_():
            self.Mix = dialog.value
            
            
    
    '''
    Control layout end
    '''    
    
        
    def OnClick(self,event): 
        if event.button == 3:
            depth = event.xdata 
            #print (depth)
            dh = self.sys.Depth[-1]/self.sys.num
            n = int(depth/dh)
            
            phase = self.Phase_diagram[n][3:53]
            P=self.Phase_diagram[n][0]
            T=self.Phase_diagram[n][2]
            #print (P,T,n)
            dialog = Mineral_Params(Pressure=P,Temperature=T,Volume=phase)
            if dialog.exec_():
                    pass
       

        
    def Layout_Control(self):
        self.BTN() 
        
        self.PLOT_Button1=QPushButton(self)
        self.PLOT_Button1.setText('Plot Update')
        self.PLOT_Button1.clicked.connect(self.plot_Vp_Vs)        
        self.BTN_Widget = QWidget(self)
        self.BTN_Widget.setLayout(self.BTN_layout)        
            
        
        
        
        
    def plot_Vp_Vs(self):
        #self.ax5.set_yticklabels([])
        #self.ax5.set_xticklabels([])
        try:
            self.figure.delaxes(self.ax5)
            self.figure.delaxes(self.ax4)
        except:
            pass
        #self.ax4.set_xlabel(" ", fontname="Times New Roman",fontsize=10)
        if self.Vp[100] == 0:
            self.showdialog(message1= "Press Calculate Vp and Vs",message2="Press Calculate Vp and Vs")
            return None

        if self.plot_seismic_type_choice==0:
            self.ax3.cla()
        
        if self.plot_seismic_type_choice==1:
            self.ax3.cla()
            self.ax3.plot(PREM_Suzan.table_depth,PREM_Suzan.table_vp, 'g--',label = 'PREM')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)
            self.ax3.plot(PREM_Suzan.table_depth,PREM_Suzan.table_vs, 'g--',label = None)
            self.ax3.plot(PREM_Suzan.table_depth,PREM_Suzan.table_vp/PREM_Suzan.table_vs, 'g--',label = None)
                

         
        if self.plot_seismic_type_choice==2:
            self.ax3.cla()
            self.ax3.plot(AK135.table_depth/1000,AK135.table_vp/1000, 'c--',label = 'AK135')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)
            self.ax3.plot(AK135.table_depth/1000,AK135.table_vs/1000, 'c--',label = None)
            self.ax3.plot(AK135.table_depth[2:]/1000,(AK135.table_vp[2:]/AK135.table_vs[2:]), 'c--',label = None)

        if self.plot_seismic_type_choice==3:
            self.ax3.cla()
            self.ax3.plot(IASP91.table_depth/1000,IASP91.table_vp/1000, 'y--',label = 'IASP91')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)
            self.ax3.plot(IASP91.table_depth/1000,IASP91.table_vs/1000, 'y--',label = None)
            self.ax3.plot(IASP91.table_depth[2:]/1000,IASP91.table_vp[2:]/IASP91.table_vs[2:], 'y--',label = None)
            
            
        if self.plot_seismic_type_choice==5:
            self.ax3.cla()
            Depth=[];Vp=[];Vs=[]
            file1 = open(self.address_input_seismic,'r')
            for line in file1:
                if line[0] == '#':
                    pass
                else:
                    line = line.split()
                    Depth.append(float(line[0]))
                    Vs.append(float(line[1]))
            Depth=np.array(Depth);Vs=np.array(Vs)
            
            self.ax3.plot(Depth,Vs, 'k',label = 'User input')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)

        
        if self.plot_seismic_type_choice==4:
            self.ax3.cla()
            self.ax3.plot(PEMC.table_depth,PEMC.table_vp, 'm--',label = 'PEMC')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)
            self.ax3.plot(PEMC.table_depth,PEMC.table_vs, 'm--',label = None)
            self.ax3.plot(PEMC.table_depth,PEMC.table_vp/PEMC.table_vs, 'm--',label = None)
            
        if self.plot_seismic_type_choice==7:
            self.ax3.cla()
            self.ax3.plot(Suzan_ref.table_depth/1000.,Suzan_ref.table_vs/1000., 'm--',label = 'NA04-ref')           
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)
        if self.plot_seismic_type_choice==8:
            self.ax3.cla()
            self.ax3.plot(Suzan_nul.table_depth/1000.,Suzan_nul.table_vs/1000., 'm--',label = 'NA04-1')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)
        if self.plot_seismic_type_choice==9:
            self.ax3.cla()
            self.ax3.plot(Suzan_one.table_depth/1000.,Suzan_one.table_vs/1000., 'm--',label = 'NA04-2')
            self.ax3.legend(loc='upper right', bbox_to_anchor=(1, 1.15),
                ncol=3, fancybox=True)            
                        
        if self.plot_type_choice != 8:
            #self.ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, borderaxespad=0.)
            self.ax2.set_xlim(10,800)
        if self.plot_seismic_type_choice != 0:
            #self.ax3.legend(loc=2,ncol=3, borderaxespad=0.)
            self.ax3.set_xlim(10,800)   
            self.ax2.set_xlim(10,800)
            
       
        if self.plot_type_choice==1:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.Vp, 'r',label = 'Calculated Vp')
            self.ax2.set_ylabel('$V_p$ (km/s)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(6.5,12)
            self.ax3.set_ylim(6.5,12)
            #self.ax4.set_ylim(6.5,12)
            #self.ax5.set_ylim(6.5,12)
            self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
            
        if self.plot_type_choice==2:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.Vs, 'g',label = 'Calculated Vs' )
            self.ax2.set_ylabel('$V_s$  (km/s)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(3.5,7.5)
            self.ax3.set_ylim(3.5,7.5)
            #self.ax4.set_ylim(3.5,7.5)
            #self.ax5.set_ylim(3.5,7.5)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
            
        if self.plot_type_choice==3:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.Vp, 'r',label = 'Calculated Vp')
            self.ax2.plot(self.sys.Depth, self.Vs, 'g',label = 'Calculated Vs')
            self.ax2.set_ylabel('$V_p$ ,$V_s$  (km/s)', color='k', fontname="Times New Roman",fontsize=10) 
            self.ax2.set_ylim(3.5,13)
            self.ax3.set_ylim(3.5,13)
            #self.ax4.set_ylim(3.5,13)
            #self.ax5.set_ylim(3.5,13)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
            
        if self.plot_type_choice==4:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.Vp/self.Vs, 'r',label = 'Vp/Vs')
            self.ax2.set_ylabel('$V_p$ /$V_s$ ', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(1,2.5)
            self.ax3.set_ylim(1,2.5)
            #self.ax4.set_ylim(1,2.5)
            #self.ax5.set_ylim(1,2.5)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
            
        if self.plot_type_choice==5:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.sys.Temperature, 'r',label = 'Temperature ')
            self.ax2.set_ylabel('Temperature (K)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(self.sys.Temperature[0],self.sys.Temperature[-1])
            self.ax3.set_ylim(self.sys.Temperature[0],self.sys.Temperature[-1]) 
            #self.ax4.set_ylim(self.sys.Temperature[0],self.sys.Temperature[-1]) 
            #self.ax5.set_ylim(self.sys.Temperature[0],self.sys.Temperature[-1])            
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
            
            
        if self.plot_type_choice==6:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.K/1e9, 'g',label = 'Bulk Modulus ')
            self.ax2.set_ylabel('Bulk Modulus (GPa)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(50,300)
            self.ax3.set_ylim(50,300)
            #self.ax4.set_ylim(50,300)
            #self.ax5.set_ylim(50,300)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)

        if self.plot_type_choice==7:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.G/1e9, 'g',label = 'Shear Modulus ')
            self.ax2.set_ylabel('Shear Modulus (GPa)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(40,220)
            self.ax3.set_ylim(40,220)
            #self.ax4.set_ylim(40,220)
            #self.ax5.set_ylim(40,220)
            self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
        
        if self.plot_type_choice==8:
            self.ax2.cla()
            self.ax2.set_ylim(3.5,13)
            self.ax3.set_ylim(3.5,13)
            #self.ax4.set_ylim(3.5,13)
            #self.ax5.set_ylim(3.5,13)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
            

        if self.plot_type_choice==9:
            self.ax2.cla()
            self.ax2.plot(self.sys.Depth, self.Rho, 'r',label = 'Density')
            self.ax2.set_ylabel('Density (kg/m$^3$)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(3000,6000)
            self.ax3.set_ylim(3000,6000)
            #self.ax4.set_ylim(3000,6000)
            #self.ax5.set_ylim(3000,6000)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)

        if self.plot_type_choice==10:
            self.ax2.cla()
            self.ax3.cla()
            self.ax2.plot(self.sys.Depth, self.dlnVp*1e2, 'r',label = 'dlnVp/dD')
            self.ax2.plot(self.sys.Depth, self.dlnVs*1e2, 'g',label = 'dlnVs/dD')            
            self.ax2.set_ylabel('dlnV/dD ($10^2$km$^-$$^1$)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(0,2)
            self.ax3.set_ylim(0,2)
            #self.ax4.set_ylim(0,2)
            #self.ax5.set_ylim(0,2)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)

        if self.plot_type_choice==11:
            self.ax2.cla()
            self.ax3.cla()
            self.ax2.plot(self.sys.Depth, self.dlnVpP*1e8, 'r',label = 'dlnVp/dP')
            self.ax2.plot(self.sys.Depth, self.dlnVsP*1e8, 'g',label = 'dlnVs/dP')            
            self.ax2.set_ylabel('dlnV/dP ($10^8$Bar$^-$$^1$)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(0,10)
            self.ax3.set_ylim(0,10)
            #self.ax4.set_ylim(0,10) 
            #self.ax5.set_ylim(0,10)           
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)

        if self.plot_type_choice==12:
            self.ax2.cla()
            self.ax3.cla()
            self.ax2.plot(self.sys.Depth, self.dlnVpT*1e6, 'r',label = 'dlnVp/dT')
            self.ax2.plot(self.sys.Depth, self.dlnVsT*1e6, 'g',label = 'dlnVs/dT')            
            self.ax2.set_ylabel('dlnV/dT ($10^6$k$^-$$^1$)', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(-3,3)
            self.ax3.set_ylim(-3,3)
            #self.ax4.set_ylim(-3,3)
            #self.ax5.set_ylim(-3,3)
            #self.ax3.set_xlim(10,800)
            self.ax2.set_xlim(10,800)

        if self.plot_type_choice==13:
            self.ax2.cla()
            self.ax3.cla()
            self.ax2.plot(self.sys.Depth,self.sys.Pressure/1e4, 'r',label = 'Pressure (GPa)')          
            self.ax2.set_ylabel('GPa', color='k', fontname="Times New Roman",fontsize=10)
            self.ax2.set_ylim(0,30)
            self.ax3.set_ylim(0,30)
            self.ax3.set_xlim(10,800)
            #self.ax4.set_xlim(10,800)
            #self.ax5.set_xlim(10,800)
            self.ax2.set_xlim(10,800)
             
        self.ax2.legend(loc='upper left', bbox_to_anchor=(-0.15,1.15),
          ncol=3, fancybox=True,)
        self.qmc.draw() 
            
    def plot_function(self):
        '''
        This plot the phase diagram over whole upper mantle.
        '''
        
        self.figure.clear()
        self.ax = self.qmc.figure.add_subplot(111) 
        self.ax2 = self.ax.twinx()
        self.ax3 = self.ax.twinx()
        self.ax4 = self.ax.twiny()
        self.ax5 = self.ax.twinx()
        

        mineral_list=['P','D','T','an','ab','sp','hc','en','fs','mgts','odi','hpcen','hpcfs','di','he','cen','cats','jd','py','al','gr','mgmj','jdmj', \
                      'capv','fo','fa','mgwa','fewa','mgri','feri','mgil','feil','co','mgpv','fepv','alpv','mppv','fppv','appv','mgcf','fecf', \
                      'nacf','pe','wu','qtz','coes','st','apbo','ky','neph','OL','WA','RI','Vp','Vs']

        name=['sp','cpx','opx','ol','plg','ak','fp','hpcpx','gt','wa','ri','pv','ppv','cf','q','coe','st','capv']
        


        
        global_Mineral_V=dict()
        
        for i,name in  enumerate(mineral_list, start=0):
            global_Mineral_V[name]=self.Phase_diagram[:,i]
        self.global_Mineral_V = global_Mineral_V

        global_m=len(global_Mineral_V['ab'])
        Mineral_V_Plot=dict()

        for i in name:   
            Mineral_V_Plot[i]=np.zeros(global_m)


            
        Mineral_V_Plot['sp']=global_Mineral_V['sp']+global_Mineral_V['hc']
        Mineral_V_Plot['cpx']=global_Mineral_V['di']+global_Mineral_V['he']+global_Mineral_V['cen']+global_Mineral_V['cats']+global_Mineral_V['jd']
        Mineral_V_Plot['opx']=global_Mineral_V['en']+global_Mineral_V['fs']+global_Mineral_V['mgts']+global_Mineral_V['odi']
#==============================================================================
#         try:
#             Mineral_V_Plot['ol']=global_Mineral_V['OL']
#             Mineral_V_Plot['wa']=global_Mineral_V['WA']
#             Mineral_V_Plot['ri']=global_Mineral_V['RI']
#         except:
#==============================================================================
        Mineral_V_Plot['ol']=global_Mineral_V['fo']+global_Mineral_V['fa']
        Mineral_V_Plot['wa']=global_Mineral_V['mgwa']+global_Mineral_V['fewa']
        Mineral_V_Plot['ri']=global_Mineral_V['mgri']+global_Mineral_V['feri']
        
        Mineral_V_Plot['plg']=global_Mineral_V['ab']+global_Mineral_V['an']
        Mineral_V_Plot['fp']=global_Mineral_V['pe']+global_Mineral_V['wu']
        
        Mineral_V_Plot['hpcpx']=global_Mineral_V['hpcen']+global_Mineral_V['hpcfs']
        Mineral_V_Plot['ak']=global_Mineral_V['co']
        
        Mineral_V_Plot['gt']=global_Mineral_V['py']+global_Mineral_V['al']+global_Mineral_V['gr']+global_Mineral_V['mgmj']+global_Mineral_V['jdmj']
        Mineral_V_Plot['pv']=global_Mineral_V['mgpv']+global_Mineral_V['fepv']+global_Mineral_V['alpv']
        Mineral_V_Plot['ppv']=global_Mineral_V['mppv']+global_Mineral_V['fppv']+global_Mineral_V['appv']+global_Mineral_V['capv']
        Mineral_V_Plot['cf']=global_Mineral_V['mgcf']+global_Mineral_V['fecf']+global_Mineral_V['nacf']
        Mineral_V_Plot['q']=global_Mineral_V['qtz']
        Mineral_V_Plot['coe']=global_Mineral_V['coes']
        Mineral_V_Plot['st']=global_Mineral_V['st']
        Mineral_V_Plot['capv']=global_Mineral_V['capv']

#==============================================================================
#         if self.Discontinuity410:
#             self.ax5.plot([self.Discontinuity410,self.Discontinuity410],[0,1])
#             self.ax5.plot([self.Discontinuity520,self.Discontinuity520],[0,1])
#             self.ax5.plot([self.Discontinuity660,self.Discontinuity660],[0,1])
#         else:
#             for loop in range(global_m):
#                 if Mineral_V_Plot['wa'][loop] >=Mineral_V_Plot['ol'][loop]:
#                     self.Discontinuity410 = self.sys.Depth[loop]
#                     a = loop
#                     break
#             for loop in range(a,global_m):
#                 if Mineral_V_Plot['ri'][loop] >=Mineral_V_Plot['wa'][loop]:
#                     self.Discontinuity520 = self.sys.Depth[loop]
#                     a = loop
#                     break            
#             for loop in range(a,global_m):
#                 if Mineral_V_Plot['ri'][loop] ==0:
#                     self.Discontinuity660 = self.sys.Depth[loop]
#                     a = loop
#                     break
#             self.ax5.plot([self.Discontinuity410,self.Discontinuity410],[0,1],'g')
#             self.ax5.plot([self.Discontinuity520,self.Discontinuity520],[0,1],'g')
#             self.ax5.plot([self.Discontinuity660,self.Discontinuity660],[0,1],'g')    
#==============================================================================
            
        self.Mineral_V_Plot = Mineral_V_Plot
        istart=0;iend=0
        for i in range(len(Mineral_V_Plot['ri'])):
            if Mineral_V_Plot['ri'][i]!=0:
                istart=i;break
        for i in range(istart,len(Mineral_V_Plot['ri'])):
            if Mineral_V_Plot['ri'][i]==0:
                iend=i;break   
        for i in range(istart,iend):
            Mineral_V_Plot['ri'][i]+=Mineral_V_Plot['fp'][i]
            Mineral_V_Plot['fp'][i]=0
            
        line1=[];line2=[];line3=[];line4=[];line5=[];line6=[];line0=[];line01=np.zeros(len(global_Mineral_V['st']));
        line11=np.zeros(len(global_Mineral_V['st']));line21=np.zeros(len(global_Mineral_V['st']))
        line31=np.zeros(len(global_Mineral_V['st']));line41=np.zeros(len(global_Mineral_V['st']))
        line51=np.zeros(len(global_Mineral_V['st']));line61=np.zeros(len(global_Mineral_V['st']))
    
        def smooth(Plot):
            #return Plot
            for i in range(2,len(Plot)-3):
                Plot[i]=(Plot[i-1]+Plot[i]+Plot[i+1])/3.
            return Plot   
            
        name1=['ol','wa','ri','fp']
        for i in range(len(name1)):
            line1.append((Mineral_V_Plot[name1[i]]))
            line11+=Mineral_V_Plot[name1[i]]
    
        name2=['pv','gt']
        for i in range(len(name2)):
            line2.append((Mineral_V_Plot[name2[i]])+line21)
            line21+=Mineral_V_Plot[name2[i]]
            
        name3=['opx','hpcpx']
        for i in range(len(name3)):
            line3.append((Mineral_V_Plot[name3[i]])+line31) 
            line31+=Mineral_V_Plot[name3[i]]
    
        name4=['cpx','plg','sp','cf']
        for i in range(len(name4)):
            line4.append((Mineral_V_Plot[name4[i]])+line41) 
            line41+=Mineral_V_Plot[name4[i]]
        
        name5=['q','coe','st','capv']
        for i in range(len(name5)):
            line5.append((Mineral_V_Plot[name5[i]])+line51) 
            line51+=Mineral_V_Plot[name5[i]]   
        

        global_D=global_Mineral_V['D']
        for i in range(len(name1)):
            if i ==1:
                pass
            else:
                self.ax.plot(global_D[:-1],line1[i][:-1]+line01[:-1],color='b')#Ol
        #ax.plot(global_D[:-10],line11[:-10],color='b')   
        for i in range(len(name2)):
            self.ax.plot(global_D[:-1],(line2[i][:-1]+line01[:-1]+line11[:-1]),color='b')
        for i in range(len(name3)):
            self.ax.plot(global_D[:-1],(line3[i][:-1]+line01[:-1]+line11[:-1]+line21[:-1]),color='b')
        for i in range(len(name4)):
            self.ax.plot(global_D[:-1],(line4[i][:-1]+line01[:-1]+line11[:-1]+line21[:-1]+line31[:-1]),color='b')
        for i in range(len(name5)):
            self.ax.plot(global_D[:-1],(line5[i][:-1]+line01[:-1]+line11[:-1]+line21[:-1]+line31[:-1]+line41[:-1]),color='b')
       
   
        a=250.
        if Mineral_V_Plot['ol'][int(50/a*global_m)]>=0.1:
            self.ax.text(global_D[int(50/a*global_m)], 0.5*Mineral_V_Plot['ol'][int(50/a*global_m)], 'ol', fontsize=12, picker=True)
        if Mineral_V_Plot['wa'][int(140/a*global_m)]>=0.1:
            self.ax.text(global_D[int(140/a*global_m)], 0.5*Mineral_V_Plot['wa'][int(140/a*global_m)], 'wa', fontsize=12, picker=True)
        if Mineral_V_Plot['ri'][int(180/a*global_m)]>=0.1:
            self.ax.text(global_D[int(180/a*global_m)], 0.5*Mineral_V_Plot['ri'][int(180/a*global_m)], 'ri', fontsize=12, picker=True)
            #ax.text(global_D[int(280/a*global_m)], 0.5*Mineral_V_Plot['fp'][int(280/a*global_m)], 'fp', fontsize=12)
        self.ax.text(global_D[int(160/a*global_m)], 0.5*Mineral_V_Plot['gt'][int(160/a*global_m)]+line11[int(160/a*global_m)], 'gt', fontsize=12, picker=True)
        #ax.text(global_D[int(280/a*global_m)], 0.5*Mineral_V_Plot['pv'][int(280/a*global_m)]+line11[int(280/a*global_m)], 'pv', fontsize=12)
        if Mineral_V_Plot['opx'][int(10/a*global_m)] >=0.05:
            self.ax.text(global_D[int(10/a*global_m)], 0.7*Mineral_V_Plot['opx'][int(40/a*global_m)]+line11[int(40/a*global_m)], 'opx', fontsize=12, picker=True)
        if Mineral_V_Plot['hpcpx'][int(130/a*global_m)] >=0.02:
            self.ax.text(global_D[int(115/a*global_m)], 0.8*Mineral_V_Plot['hpcpx'][int(130/a*global_m)]+line11[int(130/a*global_m)]+Mineral_V_Plot['gt'][int(40/a*global_m)], 'hpcpx', fontsize=12, picker=True) 
        self.ax.text(global_D[int(40/a*global_m)], 0.5*Mineral_V_Plot['cpx'][int(40/a*global_m)]+line11[int(40/a*global_m)]+Mineral_V_Plot['gt'][int(40/a*global_m)]+Mineral_V_Plot['opx'][int(40/a*global_m)], 'cpx', fontsize=12, picker=True) 
        if Mineral_V_Plot['st'][int(150/a*global_m)]>=0.02:
            self.ax.text(global_D[int(40/a*global_m)], 1-0.7*Mineral_V_Plot['st'][int(150/a*global_m)], 'co', fontsize=12, picker=True)
            self.ax.text(global_D[int(150/a*global_m)], 1-0.7*Mineral_V_Plot['st'][int(150/a*global_m)], 'st', fontsize=12, picker=True)
        if Mineral_V_Plot['capv'][int(230/a*global_m)]>=0.02:
            self.ax.text(global_D[int(230/a*global_m)], 1-0.7*Mineral_V_Plot['capv'][int(230/a*global_m)], 'capv', fontsize=12, picker=True)

        if Mineral_V_Plot['pv'][int(230/a*global_m)]>=0.02:
            self.ax.text(global_D[int(230/a*global_m)], 1-0.7*Mineral_V_Plot['pv'][int(230/a*global_m)], 'pv', fontsize=12, picker=True)
        if Mineral_V_Plot['fp'][int(230/a*global_m)]>=0.02:
            self.ax.text(global_D[int(230/a*global_m)], 0.7*Mineral_V_Plot['fp'][int(230/a*global_m)], 'fp', fontsize=12, picker=True)         
            

        self.ax4.set_xlim(GUI.sys.Pressure[4]/1e4,GUI.sys.Pressure[-1]/1e4)
        self.ax4.plot(GUI.sys.Pressure/1e4,np.ones(len(GUI.sys.Pressure)))
        self.ax4.set_xlabel("Pressure (GPa)", fontname="Times New Roman",fontsize=10)
        
        self.ax5.set_xlim(global_D[4],global_D[-1])
        self.ax5.set_ylim(0,1) 
        #self.ax5.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)       
            
        self.ax.set_xlim(global_D[4],global_D[-1])
        self.ax.set_ylim(0,1)    
        self.ax.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
        self.ax.set_ylabel('Phase proportion (volume%)', fontname="Times New Roman",fontsize=10)
        
        #self.ax4.xaxis.set_label_coords(0.5,0.95)
        self.qmc.draw()  

    
    
    
    
    def StatusBar(self):
        self.status = self.statusBar()
        self.status.setSizeGripEnabled(False)
        self.status.showMessage("Ready", 5000)        
    
        
        
    '''
    MenuBar and It's actions
    
    '''
    def createAction(self, text, slot=None, shortcut=None, icon=None,
                     tip=None, checkable=False, signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/{0}.png".format(icon)))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
            #self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action    
    
    def MenuBar(self):
        fileRestart = self.createAction("&Restart",self.filerestart,
                "filerestart","Restart")
        fileOpen = self.createAction("&Open",self.Open,
                "Open","Open")        
        fileExportAction = self.createAction("Export Phase Proportion", self.fileSaveAs,
                 "Export phase proportion", "Export phase proportion")    
        fileExportVPVSAction = self.createAction("Export Vp,Vs", self.fileExport,
                 "Export velosity", "Export velosity")  
        fileSaveAction = self.createAction("Save", self.fileSave,
                 "Save", "Save")  
        fileSaveTemperature = self.createAction("Save Temperature", self.fileSaveT,
                 "T", "T") 
        
        fileQuitAction = self.createAction("&Quit", self.fileClose,
                "Ctrl+Q", "filequit", "Close the application")

        fileview = self.createAction("&View/Modify Thermodynamics Data", self.mineral_params,
                 "View/Modify Thermodynamics Data", "View/Modify Thermodynamics Data")           
        fileEdit = self.createAction("&Mineral Plot", self.Mineralplot,
                 "Mineral plot", "Mineral plot")        
        fileTemperatureEdit = self.createAction("&Change temperature profile", self.ChangeTemperature,
                 "Change Temperature", "Change Temperature") 
        fileExportfigure = self.createAction("&Export figure", self.exportfigure,
                 "Export figure", "Export figure") 
        fileChangexaxes = self.createAction("&Change depth to pressure", self.Changexaxes,
                 "Change depth to pressure", "Change depth to pressure") 
        
        helpAboutAction = self.createAction("&About Seismic Calculator",
                self.helpAbout)
        helpHelpAction = self.createAction("&Help", self.helpHelp,
                QKeySequence.HelpContents)        
        
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(fileOpen)
        self.fileMenu.addAction(fileRestart)
        self.fileMenu.addAction(fileExportAction)
        self.fileMenu.addAction(fileExportVPVSAction)
        self.fileMenu.addAction(fileSaveAction)
        self.fileMenu.addAction(fileSaveTemperature)
        
        self.separatorAct = self.fileMenu.addSeparator()
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(fileQuitAction)
        self.menuBar().addSeparator()
        
        self.fileMenu = self.menuBar().addMenu("&Edit")
        self.fileMenu.addAction(fileview)
        self.fileMenu.addAction(fileEdit)
        self.fileMenu.addAction(fileTemperatureEdit)
        self.fileMenu.addAction(fileExportfigure)
        #self.fileMenu.addAction(fileChangexaxes)
        
        self.helpMenu = self.menuBar().addMenu("&Help")
        self.helpMenu.addAction(helpAboutAction)
        self.helpMenu.addAction(helpHelpAction)
     

    def Changexaxes(self):
        Pressure = GUI.sys.Pressure/1e4
        pass
               

        
    def okToContinue(self):
        if self.dirty:
            reply = QMessageBox.question(self,
                    "Image Changer - Unsaved Changes",
                    "Save unsaved changes?",
                    QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                return self.fileSaveAs()
        return True

    def Mineralplot(self):
        if self.Vp[100] == 0:
            self.showdialog(message1= "Please finish step 1-7",message2="Please finish step 1-7")
            return None

        if self.Check_water.isChecked() != True:
            water = False
        else:
            water = True
        if self.OLOL is None:
            self.OLOL = fo
            self.WAWA = mgwa
            self.RIRI = mgri
            
        window = GUI_tree(c2c ,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,Pv,Ring,Sp,Wad,self.OLOL,self.WAWA,self.RIRI,water)
        #print ('wf')
        if window.show():
            pass


    def Change_PT(self,pressure=None,temperature=None):
            self.progress = QProgressDialog("Running","Cancel",0,100,self) 
            self.progress.setWindowTitle('Please wait...')
            self.progress.setWindowModality(Qt.WindowModal)
            self.progress.setMinimumDuration(0)
            self.progress.setValue(0)
            self.progress.canceled.connect(self.progress.close)
            for i in range(10):
                self.progress.setValue(i)
            
            #self.sys.composition_profile(composition=cc)
            self. sys.Phase_diagram_temperature_profile(temperature)   
            #print ('calculate')
            
            for i in range(10,50):
                self.progress.setValue(i)
                
            for number,phase in enumerate(self.sys.Phase_diagram):
                self.fo[number]=phase[24]/(phase[24]+phase[25]+1e-11);
                self.mgwa[number]=phase[26]/(phase[26]+phase[27]+1e-11);
                self.mgri[number]=phase[28]/(phase[28]+phase[29]+1e-11);
            #print ('wf')
            
            for i in range(50,90):
                self.progress.setValue(i)
                
            self.Phase_diagram = np.array(self.sys.Phase_diagram)
            self.plot_function()
            for i in range(90,101):
                self.progress.setValue(i)        

    
    def ChangeTemperature(self):
        if PYQT==5:
            fileName,_ = QFileDialog.getOpenFileName(self, 'Open file', '/home')   
        else:
            fileName = QFileDialog.getOpenFileName(self, 'Open file', '/home') 
        D = [];T = []
        if fileName:
            file1 = open(fileName,'r')
            #print ('open')
            for line in file1:
                if line[0] == '#':
                    pass
                else:
                    line = line.split()
                    D.append(line[0]);T.append(line[1])
            file1.close()
            extrapolator = UnivariateSpline( D,T, k=1 )
            temperature  = extrapolator(np.linspace(0,800,300))
            self.Change_PT(temperature=temperature)

    def ChangePressureTemperature(self):
        try:
            fileName,_ = QFileDialog.getOpenFileName(self, 'Open file', '/home')   
        except:
            fileName = QFileDialog.getOpenFileName(self, 'Open file', '/home') 
        D = [];T = [];P=[]
        if fileName:
            file1 = open(fileName,'r')
            #print ('open')
            for line in file1:
                line = line.split()
                if line[0] == '#' or line[0][0]=='#':
                    pass
                else:
                    D.append(line[0]);T.append(line[1]);P.append(line[2])
            file1.close()
            extrapolator = UnivariateSpline( D,T, k=1 )
            temperature  = extrapolator(np.linspace(0,800,300))
            self.Change_PT(temperature)
                    
    def exportfigure(self):
        options = QFileDialog.Options()
        fileName = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;", options=options)   
        if PYQT == 5:
            fileName = fileName[0]
        else:
            pass
        self.figure.savefig(fileName+'.jpg',dpi=400)
        #print ('done')
        #self.fileExport()
        #self.fileSave()

        

        
    def filerestart(self):
        self.dirty = False
        self.filename = None
        self.TT = None
        self.model_choice = None
        self.model_choice_message = None
        self.sys = Phase_diagram(300) 
        self.textEdit.setText("")  
        self.average = None
        self.AverageSchemesButton.setCurrentIndex(0)
        self.BTN_Mantle_cambo.setCurrentIndex(0)
        self.BTN_NAMs_cambo.setCurrentIndex(0)        
        self.Plot_type.setCurrentIndex(0)
        self.Plot_seismice_type.setCurrentIndex(0)
        self.ax.cla()
        self.ax2.cla()
        self.ax3.cla()     
        self.qmc.draw()  
        self.Vp=np.zeros(len(self.sys.Depth))
        self.Vs=np.zeros(len(self.sys.Depth))
        self.waterusercontorl={
            'OL_function':'4.0*1e-6*D*D',
            'WA_function':'3.3',
            'RI_function':'1.7',
            'OLcoefficient':0,
            'WAcoefficient':0,
            'RIcoefficient':0,
            }
        self.AverageSchemesButton.setCurrentIndex(self.average)
        global olivine,wadsleyte,ringwoodite
        #print ('wf')
        olivine=Regression('Olivine')
        wadsleyte=Regression('Wadsleyite')
        ringwoodite=Regression('Ringwoodite') 
        #print ('wf')
        
        

    def Open(self):
#==============================================================================
#         try:
#             fileName,_  = QFileDialog.getOpenFileName(self, 'Open file', '/home')   
#         except:
#==============================================================================
        fileName = QFileDialog.getOpenFileName(self, 'Open file', '/home') 
        if PYQT == 5:
            fileName = fileName[0]
        else:
            pass
        if fileName:
            file1 = open(fileName,'r')
            for num,line in enumerate(file1):
                #print (num,line)
                if num==0:
                    line = line.split('#')
                    self.model_choice = line[0]
                    self.message = line[1]
                    self.status.showMessage(self.message)  
                    
                    self.savemessage = line[2]
                    self.textEdit.setText(line[3])
                    self.TT = float(line[3])
                    self.BTN_Mantle_cambo.setCurrentIndex(self.BTN_Mantle_cambo.findText(self.savemessage))
                    self.model_choice_message = line[4]           
                    if self.model_choice == 'Mechanical Mixture':
                        a1,a2,aa=line[5],line[6],line[7]
                        self.Type = [a1,a2,float(aa)]           
                    else:
                        self.Type=line[5]
                    #print ('1')
                if num==1:
                    line = line.split('#')
                    if line[0] == 'lalaland':
                        self.stringOL=None
                    else:
                        self.stringOL=[line[0],line[1],line[2],True]
                    self.olflag = line[4]
                    #print ('2')
                if num==2:
                    line = line.split('#')
                    if line[0] == 'lalaland':
                        self.stringWA=None
                    else:
                        self.stringWA=[line[0],line[1],line[2],True]  
                    self.waflag = line[4]
                    ##print ('3')
                if num==3:
                    line = line.split('#')
                    if line[0] == 'lalaland':
                        self.stringRI=None
                    else:
                        self.stringRI=[line[0],line[1],line[2],True]
                    self.riflag = line[4]
                    #print ('4')
                if num==4:
                    line = line.split('#')
                    self.waterusercontorl['OL_function'] = line[0]
                    self.waterusercontorl['WA_function'] = line[1]
                    self.waterusercontorl['RI_function'] = line[2]
                    self.waterusercontorl['OLcoefficient'] = float(line[3])
                    self.waterusercontorl['WAcoefficient'] = float(line[4])
                    self.waterusercontorl['RIcoefficient'] = float(line[5])
                    #print ('5')
                if num==5:
                    line= line.split('#')
                    self.average = int(line[0])
                    #print ('6')
                    
                    
                if num==6:
                    line = line.split('#')
                    self.plot_type_choice = int(line[0])
                    self.plot_seismic_type_choice = int(line[1])
                    #aa = line[3]
                    #print ('7')
                    
                    
                if num==7:
                    line = line.split('#')
                    if line[0] == 'True':
                        self.Check_water.setChecked(True)
                    else:
                        self.Check_water.setChecked(False)
                    #print ('8')
                        
            file1.close()
            #print ('wf')
            self.Select_Model()
            self.calculate(Pass=True)
            self.dialog= water_profile(Depth =self.sys.Depth,Pressure =self.sys.Pressure, Temperature = self.sys.Temperature,usercontrol=self.waterusercontorl)
            self.dialog.plot()  
            self.AverageSchemesButton.setCurrentIndex(self.average)
            self.plot_calculate_velosity()
            self.plot_Vp_Vs()
            self.Plot_seismice_type.setCurrentIndex(self.plot_seismic_type_choice+1)
            self.Plot_type.findText(self.Plot_typelist[self.plot_type_choice])            
            
            
    def fileSave(self):
        if self.Vp[100] == 0:
            self.showdialog(message1= "Please finish step 1-7",message2="Please finish step 1-7")
            return None
        
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        #try:
        fileName = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)   
        #except:
        #    fileName, _= QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)   
        if PYQT == 5:
            fileName = fileName[0]
        else:
            pass
        if fileName:
            file1 = open(fileName,'w')
            file1.write(self.model_choice);file1.write('#')
            file1.write(self.message);file1.write('#')
            file1.write(self.savemessage);file1.write('#')   
            file1.write(str(self.TT));file1.write('#')
            file1.write(self.model_choice_message);file1.write('#')
            if self.model_choice =='Mechanical Mixture':
                file1.write(self.Type[0]);file1.write('#')
                file1.write(self.Type[1]);file1.write('#')
                file1.write(str(self.Type[2]));file1.write('#')
            else:
                file1.write(self.Type);file1.write('#')    
            file1.write('\n')
            #save OL WA Ri settings
            try:
                file1.write(self.stringOL[0][:]);file1.write('#');file1.write(self.stringOL[1][:]);file1.write('#');file1.write(self.stringOL[2][:]);file1.write('#');file1.write('True');file1.write('#');
               
            except:
                file1.write('lalaland');file1.write('#');file1.write('lalaland');file1.write('#');file1.write('lalaland');file1.write('#');file1.write('False');file1.write('#');
            file1.write(self.olflag);file1.write('#');file1.write('\n')

            
            try:
                file1.write(self.stringWA[0][:]);file1.write('#');file1.write(self.stringWA[1][:]);file1.write('#');file1.write(self.stringWA[2][:]);file1.write('#');file1.write('True');file1.write('#');
            except:
                file1.write('lalaland');file1.write('#');file1.write('lalaland');file1.write('#');file1.write('lalaland');file1.write('#');file1.write('False');file1.write('#');               
            file1.write(self.waflag);file1.write('#');file1.write('\n')
            
        
            try:
                file1.write(self.stringRI[0][:]);file1.write('#');file1.write(self.stringRI[1][:]);file1.write('#');file1.write(self.stringRI[2][:]);file1.write('#');file1.write('True');file1.write('#');
            except:
                file1.write('lalaland');file1.write('#');file1.write('lalaland');file1.write('#');file1.write('lalaland');file1.write('#');file1.write('False');file1.write('#');                
            file1.write(self.riflag);file1.write('#');file1.write('\n')
            
            #save water control settings
            file1.write(self.waterusercontorl['OL_function']);file1.write('#');
            file1.write(self.waterusercontorl['WA_function']);file1.write('#');
            file1.write(self.waterusercontorl['RI_function']);file1.write('#');
            file1.write(str(self.waterusercontorl['OLcoefficient']));file1.write('#');
            file1.write(str(self.waterusercontorl['WAcoefficient']));file1.write('#');
            file1.write(str(self.waterusercontorl['RIcoefficient']));file1.write('#');
            file1.write('\n')
            
            file1.write(str(self.average));file1.write('#');
            file1.write('\n')
            
            file1.write(str(self.plot_type_choice));file1.write('#');
            file1.write(str(self.plot_seismic_type_choice));file1.write('#');
            #file1.write(str(self.Plot_typelist[self.plot_seismic_type_choice]));file1.write('#');
            file1.write('\n')
            
            if self.Check_water.isChecked() == True:
                file1.write('True');file1.write('#');
            else:
                file1.write('False');file1.write('#');
            file1.write('\n')    
            
            file1.close()                

            
    def fileSaveT(self):
        options = QFileDialog.Options()
        fileName = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)   
        if PYQT == 5:
            fileName = fileName[0]
        else:
            pass
        if fileName:
            file1 = open(fileName,'w') 
            for i in range(len(self.global_Mineral_V['fo'])):
                file1.write('{:.8}'.format(self.sys.Depth[i]));file1.write(' ')
                file1.write('{:.8}'.format(self.sys.Temperature[i]));file1.write(' ')
                file1.write('\n')                    
            file1.close()
    
    def fileSaveAs(self):
        mineral_list=['P','D','T','an','ab','sp','hc','en','fs','mgts','odi','hpcen','hpcfs','di','he','cen','cats','jd','py','al','gr','mgmj','jdmj', \
          'capv','fo','fa','mgwa','fewa','mgri','feri','mgil','feil','co','mgpv','fepv','alpv','mppv','fppv','appv','mgcf','fecf', \
          'nacf','pe','wu','qtz','coes','st','apbo','ky','neph','OL','WA','RI']  
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        #try:
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
            for name in mineral_list:
                file1.write('{:10.5}'.format(name))
                file1.write('  ')
            file1.write('\n')           
            for i in range(len(self.global_Mineral_V['fo'])):
                for name in mineral_list:
                    file1.write('{:.8}'.format(self.global_Mineral_V[name][i]))
                    file1.write('  ')
                file1.write('\n') 
            file1.close()
            #print ('save')
            
    def fileExport(self):
        list1=['P(Bar)','D(km)','T(K)','Vp(km/s)','Vs(km/s)','Density(kg/m3)']  
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        #try:
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
            for name in list1:
                file1.write(name)
                file1.write('  ')
            file1.write('\n')  
            for i in range(len(self.global_Mineral_V['fo'])):
                file1.write('{:.8}'.format(self.sys.Pressure[i]));file1.write(' ')
                file1.write('{:.8}'.format(self.sys.Depth[i]));file1.write(' ')
                file1.write('{:.8}'.format(self.sys.Temperature[i]));file1.write(' ')
                file1.write('{:.8}'.format(self.Vp[i]));file1.write(' ')
                file1.write('{:.8}'.format(self.Vs[i]));file1.write(' ')
                file1.write('{:.8}'.format(self.Rho[i]));file1.write(' ')
                file1.write('\n')                    
            file1.close()
            #print ('save')
            
    def fileClose(self):
        sys.exit()

    def helpAbout(self):
        QMessageBox.about(self, "About Seismic calculator",
                """<b>HyMaTZ</b> v {0}
                <p>Copyright &copy; .
                <p> This application is for visualizing changes in the seismic velocities, calculated using
ray theory, due to water content and iron content in the nominally anhydrous minerals.
                <p>Python {1} - Qt {2} - PyQt {3} on {4}""".format(
                __version__, platform.python_version(),
                QT_VERSION_STR, PYQT_VERSION_STR,
                platform.system()))


    def helpHelp(self):
        pass
        

    
if __name__ == "__main__":
    print (os.path.dirname(os.path.realpath(__file__)))
    #from Mineral_Physics.Solidsolution import c2c ,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,Pv,Ring,Sp,Wad
    qapp = QApplication(sys.argv)
    GUI = MainWindow()
    GUI.show()
    sys.exit(qapp.exec_())
    #Phase_diagram = np.array(GUI.BTN.sys.Phase_diagram)
