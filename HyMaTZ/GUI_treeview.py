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
except:
    from PyQt5.QtWidgets import (QVBoxLayout,QGroupBox,QTreeWidget,QTreeWidgetItem,QProgressDialog,QApplication,QAction,QActionGroup,QDockWidget,QFileDialog,QFrame,QInputDialog,QLabel,QListWidget,QMainWindow
                                 ,QCheckBox,QMessageBox,QHBoxLayout,QSpinBox,QWidget,QGridLayout,QTableView,QDialog,QPushButton,QSplitter,QLineEdit,QComboBox,QDialog)
    from PyQt5.QtGui import (QColor,QIcon,QImage,QImageReader,QImageWriter,QPainter,QPixmap,QStandardItemModel,
                             QStandardItem,QKeySequence)
    from PyQt5.QtCore import (QThread,pyqtSignal,PYQT_VERSION_STR,QFile,QFileInfo,QSettings,QT_VERSION_STR,QTimer,QVariant,Qt)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar     
from matplotlib.figure import Figure 
from Mineral_Physics.Solidsolution import OL_water,WA_water,RI_water

  
class GUI_tree(QDialog):
    def __init__(self,c2c,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,Pv,Ring,Sp,Wad,OL,WA,RI,water=False):
        super(GUI_tree, self).__init__()
        if water == False:
            self.MineralList = [c2c ,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,Pv,Ring,Sp,Wad]
        else:
            OLwater = OL_water(OL)
            WAwater = WA_water(WA)
            RIwater = RI_water(RI)
            self.MineralList = [c2c,CF,Cpx,Gt,Aki,Wus,Opx,Pl,Ppv,Pv,Sp,OLwater,WAwater,RIwater]
        self.setupUi()
        self.Tree()

  
    def Tree(self):

  
        # Effects.
        self.effectsTreeWidget.clear()
  
        for i in self.MineralList:
            parent = QTreeWidgetItem(self.effectsTreeWidget)
            parent.setText(0, i.name)
            parent.setFlags(parent.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)
            for x,_ in i.endmembers:
                child = QTreeWidgetItem(parent)
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
                child.setText(0, x.name)
                child.setText(1, x.formula)
                child.setCheckState(0, Qt.Unchecked)

  
    def setupUi(self):
        self.setupBackendBox() 
        

        self.Main = QHBoxLayout(self)
        self.Main_Splitter = QSplitter(Qt.Horizontal)
        self.Plot_Widget=QWidget()        
        self.layout_fig = QVBoxLayout(self.Plot_Widget)
        self.figure = Figure(dpi=200)      
        self.qmc = FigureCanvas(self.figure)  
        self.ntb = NavigationToolbar(self.qmc,self)      
        self.layout_fig.addWidget(self.qmc)
        self.layout_fig.addWidget(self.ntb)
  
        self.Main_Splitter.addWidget(self.backendBox)
        self.Main_Splitter.addWidget(self.Plot_Widget)
        self.Main.addWidget(self.Main_Splitter) 
        self.resize(1600,800)
        
        self.setWindowTitle("Mineral selection")
  
    def setupBackendBox(self):
  
        #self.effectsLabel = QLabel("Available Audio Effects:")
  
        headerLabels = ( "Name", "Formula")
  
        self.effectsTreeWidget = QTreeWidget()
        self.effectsTreeWidget.setHeaderLabels(headerLabels)
        self.effectsTreeWidget.setColumnCount(2)

        self.Plot_Combobox=QComboBox()
        self.PLOT=['Choose one','clean','Vp','Vs','Vp,Vs','Rho']
        self.Plot_Combobox.addItems(self.PLOT)
        self.Plot_Combobox.setCurrentIndex(self.Plot_Combobox.findText("Choose one"))
        self.Plot_Combobox.currentIndexChanged.connect(self.Plot_option)
  
        self.plot = QPushButton('Plot')
        self.plot.clicked.connect(self.MakePlot)
        
        
        
    
        layout = QVBoxLayout()
        layout.addWidget(self.effectsTreeWidget)
        layout.addWidget(self.Plot_Combobox)
        layout.addWidget(self.plot)        
        #layout.setRowStretch(3, 100)
  
        self.backendBox = QGroupBox()
        self.backendBox.setLayout(layout)

    def Plot_option(self,a):
        self.plotoption = a        


    def MakePlot(self):
        a = self.find_checked()
        self.ax = self.figure.add_subplot(111) 
        #print (a)
        if self.plotoption == 1:
            self.ax.cla()
        if self.plotoption == 2:#Vp
            self.ax.cla()  
            for num,i in enumerate(self.MineralList):
                for j in range(len(i.endmembers)):
                    if a[num][j] == True:
                        if len(i.endmembers[j][0].Depth) != 0:
                            self.ax.plot(i.endmembers[j][0].Depth,i.endmembers[j][0].Vp,'r')
                            number = int(0.5*len(i.endmembers[j][0].Depth))
                            position_x = i.endmembers[j][0].Depth[number]
                            position_y = i.endmembers[j][0].Vp[number]
                            self.ax.text(position_x,position_y,i.endmembers[j][0].name, fontsize=15)
            self.ax.set_ylabel('$V_p$ (km/s)', color='k', fontname="Times New Roman",fontsize=10)
            
        if self.plotoption == 3:#Vs
            self.ax.cla()  
            for num,i in enumerate(self.MineralList):
                for j in range(len(i.endmembers)):
                    if a[num][j] == True:
                        if len(i.endmembers[j][0].Depth) != 0:                        
                            self.ax.plot(i.endmembers[j][0].Depth,i.endmembers[j][0].Vs,'g')
                            number = int(0.5*len(i.endmembers[j][0].Depth))
                            position_x = i.endmembers[j][0].Depth[number]
                            position_y = i.endmembers[j][0].Vs[number]
                            self.ax.text(position_x,position_y,i.endmembers[j][0].name, fontsize=15)       
            self.ax.set_ylabel('$V_s$ (km/s)', color='k', fontname="Times New Roman",fontsize=10)

                        
        if self.plotoption == 4:#Vp,Vs
            self.ax.cla()   
            for num,i in enumerate(self.MineralList):
                for j in range(len(i.endmembers)):
                    if a[num][j] == True:
                        if len(i.endmembers[j][0].Depth) != 0:                        
                            self.ax.plot(i.endmembers[j][0].Depth,i.endmembers[j][0].Vp,'r')
                            self.ax.plot(i.endmembers[j][0].Depth,i.endmembers[j][0].Vs,'g')
                            number = int(0.5*len(i.endmembers[j][0].Depth))
                            position_x = i.endmembers[j][0].Depth[number]
                            position_y = i.endmembers[j][0].Vp[number]
                            self.ax.text(position_x,position_y,i.endmembers[j][0].name, fontsize=15)
                            position_y = i.endmembers[j][0].Vs[number]
                            self.ax.text(position_x,position_y,i.endmembers[j][0].name, fontsize=15)
            self.ax.set_ylabel('$V_p$,$V_s$ (km/s) ', color='k', fontname="Times New Roman",fontsize=10)
                        
        if self.plotoption == 5:#Rho
            self.ax.cla() 
            for num,i in enumerate(self.MineralList):
                for j in range(len(i.endmembers)):
                    if a[num][j] == True:
                        if len(i.endmembers[j][0].Depth) != 0:                        
                            self.ax.plot(i.endmembers[j][0].Depth,i.endmembers[j][0].Rho,'b')
                            number = int(0.5*len(i.endmembers[j][0].Depth))
                            position_x = i.endmembers[j][0].Depth[number]
                            position_y = i.endmembers[j][0].Rho[number]
                            self.ax.text(position_x,position_y,i.endmembers[j][0].name, fontsize=15)   
            self.ax.set_ylabel('Density (kg/m$^3$)', color='k', fontname="Times New Roman",fontsize=10)
        self.ax.set_xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)                
        self.qmc.draw()
        #print ('pass')
        pass

    def find_checked(self):
        checked = []
        root = self.effectsTreeWidget.invisibleRootItem()
        signal_count = root.childCount()
    
        for i in range(signal_count):
            signal = root.child(i)
            checked_sweeps = list()
            num_children = signal.childCount()
    
            for n in range(num_children):
                child = signal.child(n)
    
                if child.checkState(0) == Qt.Checked:
                    checked_sweeps.append(True)
                else:
                    checked_sweeps.append(False)
    
            checked.append( checked_sweeps)
    
        return checked
  
  
if __name__ == '__main__':
    from Mineral_Physics.Stix2011data import OL,WA,RI
    from Mineral_Physics.Solidsolution import c2c ,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,Pv,Ring,Sp,Wad,OLwater,WAwater,RIwater,OL_water,WA_water,RI_water
    app = QApplication(sys.argv)
    app.setApplicationName("Phonon Capabilities Example")
  
    window = GUI_tree(c2c ,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,Ppv,Pv,Ring,Sp,Wad,OL,WA,RI,water=True)
    window.show()
  
    sys.exit(app.exec_())