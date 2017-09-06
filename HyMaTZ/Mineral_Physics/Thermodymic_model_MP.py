# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 00:04:44 2016

@author: Fei

Update 10/3/2016
"""


import os
import numpy as np

'''
This file store Thermodynamics Model from Perplex
'''


class Harzburgite():
    '''
    This is the base class for all Thermodyanmics Models.
    '''
    def __init__(self,composition=None,name=None,n=60000,index_P=300,index_T=200,address=None):
        self.index_P=index_P
        self.index_T=index_T
        self.name=name
        self.directory = address
        self.Name=['Vs','Vp','Density','Pressure','Temperature','Cp','V','alpha','S','H']
        self.present=1
        if composition is None:
            composition =[0,0,0,0,0,0]
        self.composition={
                'CAO':composition[0],
                'AL2O3':composition[1],
                'NA2O':composition[2],
                'MGO':composition[3],
                'FEO':composition[4],
                'SIO2':composition[5],
                            }
        self.n=n
        self.dictionary={'Vs':np.zeros(self.n),
            'Vp':np.zeros(self.n),
            'Density':np.zeros(self.n),
            'Pressure':np.zeros(self.n),
            'Temperature':np.zeros(self.n),
            'Cp':np.zeros(self.n),
            'V':np.zeros(self.n),
            'alpha':np.zeros(self.n),
            'S':np.zeros(self.n),
            'H':np.zeros(self.n),
                        }
                        
                        
        try:
            self.Read_data()
        except:
            str1="No data file as "+str(self.address)
            raise NotImplementedError(str1)
        
    
    
    def Read_data(self):
        '''
        Read Model from a txt file.
        txt file contain thermodyancmis properties at difference pressure and temperature
        at each pressure and temperature txt file contain 
        Vs: S wave velocity
        Vp: P wave velocity
        Density: Density
        OlVs,WaVs,RiVs: volume fraction of Ol,WA and RI
        fo,fa,mgwa,fewa,mgri,feri: volume fraction of each minerals
        Cp,alpha and V: Heat capacity at pressurel; thermal expansion and Volume
        All units are SI unites 
        '''
        try:
            #directory =os.path.dirname(os.path.realpath('__file__'))
            self.address = os.path.join( self.directory , 'Models',self.name,self.name+'_2.txt')
            file=open(self.address)
            
        except:
            #pass
            self.address = os.path.join( 'E:\\', 'Models',self.name,self.name+'_2.txt')
            file=open(self.address)

        try:
            file.next();linecount=0
        except:
            
            next(file);linecount=0
        
        for line in file:
            line=line.split()
            #self.dictionary['Vs'][linecount]=line[10]
            self.dictionary['Pressure'][linecount]=line[0]
            self.dictionary['Temperature'][linecount]=line[1]
            self.dictionary['H'][linecount]=line[3]
            self.dictionary['Density'][linecount]=line[7]
            self.dictionary['Cp'][linecount]=line[5] 
            self.dictionary['alpha'][linecount]=line[4]
            self.dictionary['V'][linecount]=line[6] 
            self.dictionary['S'][linecount]=line[2]
            linecount+=1
        file.close()
 
    def find_pressure_index(self,pressure):#Bar
        '''
        Find index (lower bound)
        '''
        
        p= int(pressure/(300000/self.index_P))
        if p>=self.index_P:
            return self.index_P
        else:
            return p
        
        
    def find_temperature_index(self,temperature):#K
        '''
        Find index (lower bound)
        '''

        t = int((temperature-773)/(2000/self.index_T))
        if t >=self.index_T:
            return self.index_T
        else:
            return t   

       
    def return_dataindex(self,P,T):
        pressure_index=self.find_pressure_index(P)  
        temperature_index=self.find_temperature_index(T) 
        n=(int(temperature_index))*self.index_P+int(pressure_index)
        return n
        
    def __add__(self,other):
        dictionary=dict()
        for i in self.Name:
            dictionary[i]=(self.dictionary[i]+other.dictionary[i])
        Model=Harzburgite( [0,0,0,0,0,0],'Harzburgite10',240000,600,400,address=self.directory) 
        Model.dictionary=dictionary
        return Model
    
    def __mul__(self,other):
        dictionary=dict()
        self.present*=other
        for i in self.Name:
            dictionary[i]=self.dictionary[i]*other
        Model=Harzburgite( [0,0,0,0,0,0],'Harzburgite10',240000,600,400,address=self.directory) 
        Model.dictionary=dictionary            
        return Model     
    



if __name__ == "__main__":
    index_P=600;index_T=400
    n = index_P*index_T          
    Harzburgite100=Harzburgite([0.81,0.53,0.,56.51,6.07,36.07],'Harzburgite100',n=n,index_P=index_P,index_T=index_T)     
    Harzburgite10=Harzburgite( [ 13.88,  10.19 ,  2.18,  14.94  , 7.06,  51.75],'Harzburgite10',n=n,index_P=index_P,index_T=index_T)    
    #Mixture=Harzburgite( [0,0,0,0,0,0],'piclogite',n=n,index_P=index_P,index_T=index_T) 
    Mixture1=Harzburgite100*0.82#+Harzburgite100*0.12
    Mixture2=Harzburgite100*0.82+Harzburgite10*0.12
    
    from matplotlib.pylab import plt
    plt.plot(Mixture1.dictionary['V'],'g')
    plt.plot(Mixture2.dictionary['V'],'b')    

