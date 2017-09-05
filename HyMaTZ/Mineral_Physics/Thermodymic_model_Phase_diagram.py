# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 16:26:26 2016

@author: Fei
"""


import os
import numpy as np

'''
This file store Thermodynamics Model from Perplex
'''


class Harzburgite_Phase():
    '''
    This is the base class for all Thermodyanmics Models.
    '''
    def __init__(self,composition=None,name=None,n=240000,index_P=600,index_T=400,address=None):
        self.name=name
        self.present=1
        self.index_P=index_P
        self.index_T=index_T
        self.directory = address
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
        #n=60000
        self.storage=[]
                        
                        
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
            #self.directory =os.path.dirname(os.path.realpath('__Main__'))
            self.address = os.path.join( self.directory , 'Models',self.name,self.name+'_3.txt')
            file=open(self.address)
            
        except:
            #pass
            self.address = os.path.join( 'E:\\', 'Models',self.name,self.name+'_3.txt')
            file=open(self.address)

        try:
            file.next();linecount=0
        except:
            next(file);linecount=0
        
        
        for line in file:
            line=line.split()
            a=[float(i) for i in line]
            self.storage.append(a)
        self.storage = np.array(self.storage)


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
    
    def change_OLVs(self):
        Pi=0;Ti=2;
        an=3;ab=4;sp=5;hc=6;
        en=7;fs=8;mgts=9;odi=10;
        hpcen=11;hpcfs=12;
        di=13;he=14;cen=15;cats=16;
        jd=17;py=18;al=19;gr=20; mgmj=21;jdmj=22;capv=23;
        fo=24;fa=25;mgwa=26;fewa=27;mgri=28;feri=29;#OL=50;WA=51;RI=52
        mgil=30;feil=31;co=32;mgpv=33;fepv=34;alpv=35;mppv=36;fppv=37;appv=38;mgcf=39;fecf=40;
        nacf=41;pe=42;wu=43;qtz=44;coes=45;st=46;apbo=47;ky=48;neph=49; #apbo neph did not found
        
    
    def __add__(self,other):
        storage=(self.storage+other.storage)
        Model=Harzburgite_Phase( [0,0,0,0,0,0],'Harzburgite100') 
        Model.storage=storage
        return Model
    
    def __mul__(self,other):
        self.present*=other
        storage=self.storage*other
        Model=Harzburgite_Phase( [0,0,0,0,0,0],'Harzburgite100') 
        Model.storage=storage          
        return Model 









if __name__ == "__main__":
    Harzburgite10_Phase=Harzburgite_Phase( [3.424,2.462,0.436,48.196,6.268,39.206],'Pyrolite')      
    print (Harzburgite10_Phase.address)     
    #Harzburgite100_Phase=Harzburgite_Phase([0.81,0.53,0.,56.51,6.07,36.07],'Harzburgite100')     
    #Mixture_Phase = Harzburgite100_Phase*0.82 +Harzburgite10_Phase*0.18




