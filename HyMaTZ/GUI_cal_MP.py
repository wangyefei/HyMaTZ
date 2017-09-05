# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 10:18:36 2016

@author: Fei
"""
from __future__ import print_function
import os
import numpy as np
from scipy.interpolate import interp1d
import scipy as scipy

import matplotlib.pyplot as plt

from Mineral_Physics.Thermodymic_model_MP import Harzburgite
from Mineral_Physics.Thermodymic_model_Phase_diagram import Harzburgite_Phase
from scipy.optimize import minimize
from Mineral_Physics.tools import (Regression_Plane1)
from Seismic.Seismic import prem_Suzan0

from Mineral_Physics.Velocity_calculator import Velocity_calculator
from Mineral_Physics.regression import olivine,wadsleyte,ringwoodite
from Mineral_Physics.Solidsolution import (c2c,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,ppv,Pv,Ring,Sp,Wad)
from Mineral_Physics.Stix2011data import (OL,WA,RI,ab,an,sp,hc,fo,fa,mgwa,fewa,mgri,feri,en,fs,mgts,odi,di\
                          ,he,cen,cats,jd,hpcen,hpcfs,mgpv,fepv\
                          ,alpv,capv,mgil,feil,co,py ,al,gr,mgmj ,jdmj,qtz,coes,st,seif,mppv,fppv,appv,pe,wu,mgcf \
                          ,fecf ,nacf,ky,neph,OL_,WA_,RI_ )

g= scipy.constants.g
#==============================================================================
# directory =os.path.dirname(os.path.abspath('__file__'))
# address = os.path.join( directory,'Mineral_Physics','fort.56.txt')
# file = open(address)
# DD=[]
# PP=[]
# VVPP=[]
# VVSS=[]
# RRho=[]
# TT=[]
# #next(file)
# for line in file:
#     line = line.split()
#     DD.append(float(line[1]))
#     PP.append(float(line[0]))
#     VVSS.append(float(line[5]))
#     VVPP.append(float(line[6]))     
#     RRho.append(float(line[3]))  
#     TT.append(float(line[2]))
# file.close()  
#==============================================================================
      
class Phase_diagram():
    
    def __init__(self,composition=None,num=300,Model =None,Model_composition = None,index_P=600,index_T=400):
        '''
        composition should be two list 
        composition[0] == Depth
        composition[1] == composition parameters, if some depth not given, it will be consider as pyrolite
        '''
        if Model is not None:
            self.Model_Phase = Harzburgite_Phase(Model_composition,Model,address = os.path.dirname(os.path.realpath(__file__)))
            self.Model = Harzburgite(Model_composition,Model,n=index_P*index_T,index_P=index_P,index_T=index_T,address = os.path.dirname(os.path.realpath(__file__)))
        else:
            self.Model = None
            self.Model_Phase=None
                
        self.num=num
        self.Depth = np.linspace(0, 800, self.num, endpoint=True) 
        self.composition_profile(composition=None)
        self.Pressure = np.zeros(self.num)
        self.Temperature = np.zeros(self.num)
        self.Water = np.zeros(self.num)
        self.Pressure[0] = 1

    def change_model(self,Model='Harzburgite82',Model_composition=[2.93553,2.21491,0.109863,49.8509,6.16912,38.6950],index_P=600,index_T=400):
        self.Model_Phase = Harzburgite_Phase(Model_composition,Model,address = os.path.dirname(os.path.realpath(__file__)))
        self.Model = Harzburgite(Model_composition,Model,n=index_P*index_T,index_P=index_P,index_T=index_T,address = os.path.dirname(os.path.realpath(__file__)))     
            

    def change_model1(self,Model1,Model2,Model_composition1,index_P=600,index_T=400):
        Model1_Phase=Harzburgite_Phase([0,0,0,0,0,0],Model1,address = os.path.dirname(os.path.realpath(__file__)))
        Model2_Phase=Harzburgite_Phase([0,0,0,0,0,0],Model2,address = os.path.dirname(os.path.realpath(__file__)))
        self.Model_Phase = Model1_Phase*Model_composition1 +Model2_Phase*(1-Model_composition1)
        #print ('change1')
        Model1 = Harzburgite([0,0,0,0,0,0],Model1,n=index_P*index_T,index_P=index_P,index_T=index_T,address = os.path.dirname(os.path.realpath(__file__)))
        Model2 = Harzburgite([0,0,0,0,0,0],Model2,n=index_P*index_T,index_P=index_P,index_T=index_T,address = os.path.dirname(os.path.realpath(__file__)))
        #print ('change2')
        self.Model = Model1*Model_composition1 +Model2*(1-Model_composition1)

     
    def Import_layer(self,layers):
        self.Layers=layers
        
    def change_layer(self,depth):
        for i in range(len(self.Layers.layer_depth)):
            if self.Layers.layer_depth[i] >=depth:
                aa=i;break
        #print ('change')
        self.Model = self.Layers.layer[aa]
        self.Model_Phase = self.Layers.layer_phase[aa]            
#==============================================================================
#         if self.Model.name != self.Layers.layer[aa].name:
#             self.Model = self.Layers.layer[aa]
#             self.Model_Phase = self.Layers.ayer_phase[aa]
#         else:
#             pass
#==============================================================================
        
    
    def change_pressure(self,Pressure):
        self.Pressure = Pressure
    def change_temperature(self,temperature):
        self.Temperature = temperature
      
    '''
    This is the thermodyanmics properties calculations unit:bar, K
    '''
    
    def Phase_composition(self,Pressure,Temperature,Composition=1.01):
        '''
        Given Pressure and Temperature return information
        '''
        Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(Pressure,Temperature,Composition)
        Phase_composition = self.find_phase_diagram_at_PTC_condition(n,Pressure,Temperature,Composition)
        #print (Vs,Vp,Cp,alpha,Density,Volume,n)
        return Phase_composition
    

    def isoter(self,S=1000,P=10000,Composition=1.01):
        '''
        Given entropy return Temperature that most close to this value
        '''
        T_low=773;T_high=2773;T=0.5*(T_low+T_high)
        for i in range(200):
            Vs,Vp,Cp,alpha,Density,Volume,n,S1 = self.find_properties_at_PT_condition(P,T,Composition)
            if abs(S1-S)<=1:
                return T
            elif S1<S:
                T_low=T; T =0.5*(T_low+T_high)
            else:
                T_high = T; T =0.5*(T_low+T_high)
        return T
                
        
        
                
    def dP_dT(self,P,T,Composition=1.01,dp=1e2):
        '''
        Given delta pressure, return delta temperature (Adiabatic)
        dp : Bar
        dt : K
        '''
        dt_low = -50;dt_high = 50;dt = 0.5*(dt_low+dt_high);
        Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(P,T,Composition)
        for i in range(20):            
            Vs,Vp,Cp,alpha,Density,Volume,n,S1 = self.find_properties_at_PT_condition(P+dp,T+dt,Composition)
            if abs(S-S1)<=0.01:
                return dt 
            elif S1<S:
                dt_low=dt; dt =0.5*(dt_low+dt_high)
            else:
                dt_high = dt; dt =0.5*(dt_low+dt_high)
        return dt

            
    def Gibbs(self,Pressure,Temperature,Composition=1.01):
        '''
        Given Pressure and Temperature return Gibbs energy Unit KJ
        '''
        Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(Pressure,Temperature,Composition)
        Pressure *= 1e5
        Phase_composition = self.find_phase_diagram_at_PTC_condition(n,Pressure,Temperature,Composition)
        Sum=sum(Phase_composition[3:50])
        Phase = [Pl,Sp,Opx,c2c,Cpx,Gt,capv,O,Wad,Ring,Aki,Pv,ppv,CF,Wus,qtz,coes,st,ky,ky,neph]

        a=3;Gibbs=0
        for i in Phase:
            if sum(Phase_composition[a:a+len(i.endmembers)]) <=1e-4:
                pass
            else:
                Gibbs += sum(i.Gibbs_energy(P=Pressure,T=Temperature,mole_fractions=Phase_composition[a:a+len(i.endmembers)])*Phase_composition[a:a+len(i.endmembers)])
            a+=len(i.endmembers)  
        return Gibbs*1000/Sum
        
   
    def composition_profile(self,composition=None):
        '''
        this return composition profile within the whole upper mantle till 700km
        if user did not input model composition, it will assume the composition is Basalt
        '''
        num=self.num
        if composition is None:
            self.Composition = np.zeros(num)
            
        elif self.Depth[0]>=composition[0][0] and self.Depth[-1]<=composition[0][-1]:
            f2=interp1d(composition[0],composition[1])
            for i in range(num):
                    self.Composition[i] = f2(self.Depth[i])             
        else:
            #print ('wfwfwf')
            f2=interp1d(composition[0],composition[1])
                
                
            self.Composition = np.zeros(num)
            for number, i in enumerate(self.Depth):
                if i>composition[0][0]:
                    a1 =  number
                    break
                else:
                    a1=0
                    
            for number, i in enumerate(self.Depth):
                if i>composition[0][-1]:
                    a2=number-1
                    break
                else:
                    a2=number
                
            for i in range(num):
                if i>=a1 and i<=a2:
                    self.Composition[i] = f2(self.Depth[i])

          
        
    def Phase_diagram_first_pinciple(self,foot_temperature=1600,num=300,OLwater=np.zeros(300),WAwater=np.zeros(300),RIwater=np.zeros(300),selfcal=False,layers=False):#tested 
    
        self.Temperature[0] = foot_temperature
        dh = self.Depth[1]-self.Depth[0]
        self.Phase_diagram=[]
        if layers==False:
            pass
        else:
            self.change_layer(self.Depth[0])
        Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[0],self.Temperature[0],self.Composition[0],OLwater[0],WAwater[0],RIwater[0],selfcal)
        phase = self.find_phase_diagram_at_PTC_condition(n,self.Pressure[0],self.Temperature[0],self.Composition[0]) 
        phase[1]=0
        self.Phase_diagram.append(phase)
        for i in range(1,len(self.Depth)):
            if layers==False:
                pass
            else:
                self.change_layer(self.Depth[i])
            g=prem_Suzan0.gravity(self.Depth[i])
            Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[i-1],self.Temperature[i-1],self.Composition[i-1],OLwater[i-1],WAwater[i-1],RIwater[i-1],selfcal)
            dh = self.Depth[i]-self.Depth[i-1]
            dp = Density*dh*g/100
            self.Pressure[i] = self.Pressure[i-1]+ dp
            dp_dT = self.dP_dT(self.Pressure[i-1],self.Temperature[i-1],self.Composition[i-1],dp)
            self.Temperature[i] = self.Temperature[i-1]+dp_dT           
            for j in range(3):
                g=prem_Suzan0.gravity(self.Depth[i])
                Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[i-1],self.Temperature[i-1],self.Composition[i-1],OLwater[i],WAwater[i],RIwater[i],selfcal)
                Vs,Vp,Cp,alpha,Density1,Volume,n1,S = self.find_properties_at_PT_condition(self.Pressure[i],self.Temperature[i],self.Composition[i],OLwater[i],WAwater[i],RIwater[i],selfcal)
                D= (Density1+Density1)/2
                dh = self.Depth[i]-self.Depth[i-1]
                dp = D*dh*g/100
                self.Pressure[i] = self.Pressure[i-1]+ dp
                dp_dT = self.dP_dT(self.Pressure[i-1],self.Temperature[i-1],self.Composition[i-1],dp)
                self.Temperature[i] = self.Temperature[i-1]+dp_dT  
            
            phase = self.find_phase_diagram_at_PTC_condition(n,self.Pressure[i],self.Temperature[i],self.Composition[i]) 
            if phase[19] != 0:
                phase[22] += 1e-3
            phase[1] = self.Depth[i]
            phase[50] = phase[24] + phase[25]
            phase[51] = phase[26] + phase[27]
            phase[52] = phase[28] + phase[29]
            self.Phase_diagram.append(phase)
        return self.Temperature


    def Phase_diagram_temperature_profile(self,temperature=np.linspace(800,2000,300),num=300,OLwater=np.zeros(300),WAwater=np.zeros(300),RIwater=np.zeros(300),selfcal=False,layers=False):#tested
        self.Temperature = temperature
        dh = self.Depth[1]-self.Depth[0]
        self.Phase_diagram=[]
        if layers==False:
            pass
        else:
            self.change_layer(self.Depth[0])        
        Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[0],temperature[0],self.Composition[0],OLwater[0],WAwater[0],RIwater[0],selfcal)
        phase = self.find_phase_diagram_at_PTC_condition(n,self.Pressure[0],temperature[0],self.Composition[0]) 
        phase[1]=0
        self.Phase_diagram.append(phase)
        for i in range(1,len(self.Depth)):
            if layers==False:
                pass
            else:
                self.change_layer(self.Depth[i])           
            g=prem_Suzan0.gravity(self.Depth[i])
            Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[i-1],temperature[i-1],self.Composition[i-1],OLwater[i-1],WAwater[i-1],RIwater[i-1],selfcal)
            dp = Density*dh*g/100
            self.Pressure[i] = self.Pressure[i-1]+ dp       
            for j in range(3):
                g=prem_Suzan0.gravity(self.Depth[i])-0.05
                Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[i-1],temperature[i-1],self.Composition[i-1],OLwater[i],WAwater[i],RIwater[i],selfcal)
                Vs,Vp,Cp,alpha,Density1,Volume,n1,S = self.find_properties_at_PT_condition(self.Pressure[i],temperature[i],self.Composition[i],OLwater[i],WAwater[i],RIwater[i],selfcal)
                D= (Density1+Density1)/2
                dp = D*dh*g/100
                self.Pressure[i] = self.Pressure[i-1]+ dp            
            phase = self.find_phase_diagram_at_PTC_condition(n,self.Pressure[i],temperature[i],self.Composition[i]) 
            if phase[19] != 0:
                phase[22] += 1e-3
            phase[1] = self.Depth[i]
            phase[50] = phase[24] + phase[25]
            phase[51] = phase[26] + phase[27]
            phase[52] = phase[28] + phase[29]
            self.Phase_diagram.append(phase)
        return None


    def Phase_diagram_pressure_temperature_profile(self,pressure = np.linspace(1,3e5,300), temperature=np.linspace(800,2000,300),num=300,OLwater=np.zeros(300),WAwater=np.zeros(300),RIwater=np.zeros(300),selfcal=False):#tested
        self.Temperature = temperature
        dh = self.Depth[1]-self.Depth[0]
        self.Phase_diagram=[]
        phase = self.find_phase_diagram_at_PTC_condition(n,pressure[0],temperature[0],self.Composition[0]) 
        phase[1]=0
        self.Phase_diagram.append(phase)
        for i in range(1,len(self.Depth)): 
            phase = self.find_phase_diagram_at_PTC_condition(n,pressure[i],temperature[i],self.Composition[i]) 
            if phase[19] != 0:
                phase[22] += 1e-3
            phase[1] = self.Depth[i]
            phase[50] = phase[24] + phase[25]
            phase[51] = phase[26] + phase[27]
            phase[52] = phase[28] + phase[29]
            self.Phase_diagram.append(phase)
        return None

            
    def Phase_diagram_input(self):
        self.Phase_diagram=[]
        for i in range(len(self.Depth)):
            Vs,Vp,Cp,alpha,Density,Volume,n,S = self.find_properties_at_PT_condition(self.Pressure[i],self.Temperature[i],self.Composition[i])
            phase = self.find_phase_diagram_at_PTC_condition(n=n,P=self.Pressure[i],T=self.Temperature[i],composition=self.Composition[i]) 
            phase[1] = self.Depth[i]
            #phase = self.Phase_diagram_volume(phase)
            self.Phase_diagram.append(phase)            


            
    def Phase_diagram_volume(self,phase):
        
        phase[3]*=an.volume();phase[4]*=ab.volume();phase[5]*=sp.volume();phase[6]*=hc.volume()
        phase[7]*=en.volume();phase[8]*=fs.volume();phase[9]*=mgts.volume();phase[10]*=odi.volume();
        phase[11]*=hpcen.volume() ;phase[12]*=hpcfs.volume()
        phase[13]*=di.volume();phase[14]*=he.volume();phase[15]*=cen.volume();phase[16]*=cats.volume()
        phase[17]*=jd.volume();phase[18]*=py.volume();phase[19]*=al.volume();phase[20]*=gr.volume()
        phase[21]*=mgmj.volume();phase[22]*=jdmj.volume();phase[23]*=capv.volume()
        phase[24]*=fo.volume();phase[25]*=fa.volume()  
        phase[26]*=mgwa.volume();phase[27]*=fewa.volume()
        phase[28]*=mgri.volume();phase[29]*=feri.volume()
        phase[30]*=mgil.volume();phase[31]*=feil.volume();phase[32]*=co.volume()
        phase[33]*=mgpv.volume();phase[34]*=fepv.volume();phase[35]*=alpv.volume()
        phase[36]*=mppv.volume();phase[37]*=fppv.volume();phase[38]*=appv.volume()
        phase[39]*=mgcf.volume();phase[40]*=fecf.volume();phase[41]*=nacf.volume()
        phase[42]*=pe.volume();phase[43]*=wu.volume()
        phase[44]*=qtz.volume();phase[45]*=coes.volume()
        phase[46]*=st.volume();phase[47]*=0#apbo.volume()
        phase[48]*=ky.volume();phase[49]*=neph.volume()
        Sum=sum(phase[3:50])
        for i in range(3,50):
            phase[i]/=Sum
            #print (phase[i])
        
        phase[50]=phase[24]+phase[25];
        phase[51]=phase[26]+phase[27];
        phase[52]=phase[28]+phase[29];
        return phase
        
   
        
    

    
    def find_phase_diagram_at_PTC_condition(self,n,P=None,T=None,composition=None):
        n=int(n)
        return np.array(self.Model_Phase.storage[n])

        

    def find_properties_at_PT_condition(self,pressure,temperature,composition,OLwater=0,WAwater=0,RIwater=0,selfcal=False):
        try:
            P_index = self.Model.index_P #300
        except:
            P_index = 600
        pressure_index=self.find_pressure_index(pressure)  
        temperature_index=self.find_temperature_index(temperature) 
        n=(int(temperature_index))*P_index +int(pressure_index)   
        
        Vs=self.find_data_in_PT_MM('Vs',n,composition,pressure,temperature)
        Vp=self.find_data_in_PT_MM('Vp',n,composition,pressure,temperature)
        
        Cp=self.find_data_in_PT_MM('Cp',n,composition,pressure,temperature)
        alpha=self.find_data_in_PT_MM('alpha',n,composition,pressure,temperature)
        Density=self.find_data_in_PT_MM('Density',n,composition,pressure,temperature)
        
        if selfcal==True:
            phase = self.find_phase_diagram_at_PTC_condition(n,pressure,temperature,composition)
            Density -=  fo.Return_Rho(pressure*1e4,temperature)*phase[24]
            Density -=  fa.Return_Rho(pressure*1e4,temperature)*phase[25]
            Density -=  mgwa.Return_Rho(pressure*1e4,temperature)*phase[26]
            Density -=  fewa.Return_Rho(pressure*1e4,temperature)*phase[27]        
            Density -=  mgri.Return_Rho(pressure*1e4,temperature)*phase[28]
            Density -=  feri.Return_Rho(pressure*1e4,temperature)*phase[29]
            
            OL.Set_Water_Iron_Condition(OLwater,1-phase[24]/(phase[24]+phase[25]+1e-11))
            WA.Set_Water_Iron_Condition(WAwater,1-phase[26]/(phase[26]+phase[27]+1e-11))
            RI.Set_Water_Iron_Condition(RIwater,1-phase[28]/(phase[28]+phase[29]+1e-11))        
            Density += OL.Return_Rho(pressure*1e4,temperature)*(phase[24]+phase[25])
            Density += WA.Return_Rho(pressure*1e4,temperature)*(phase[26]+phase[27])
            Density += RI.Return_Rho(pressure*1e4,temperature) *(phase[28]+phase[29])
        
        Volume = self.find_data_in_PT_MM('V',n,composition,pressure,temperature)  
        S = self.find_data_in_PT_MM('S',n,composition,pressure,temperature)  
        return Vs,Vp,Cp,alpha,Density,Volume,n,S
            
    
            
    def find_data_in_PT_MM(self,name,n,composition,pressure,temperature):
        '''
        This function did a linear regression to find out true properties of the Rock at certain pressure and temperature
        0: lower bound of pressure and temperautre
        1: lower bound of pressure and higher bound of temperautre
        2: lower bound of temperature and higher bound of pressure 
        '''            
        try:
            P_index = self.Model.index_P #300
        except:
            P_index = 600
        z0=self.Model.dictionary[name][n]
        z1=self.Model.dictionary[name][n+P_index]
        z2=self.Model.dictionary[name][n+1]
        x0=self.Model.dictionary['Pressure'][n]
        x1=self.Model.dictionary['Pressure'][n+P_index]
        x2=self.Model.dictionary['Pressure'][n+1]
        y0=self.Model.dictionary['Temperature'][n]
        y1=self.Model.dictionary['Temperature'][n+P_index]
        y2=self.Model.dictionary['Temperature'][n+1]    

        X=[x0,x1,x2];Y=[y0,y1,y2];Z=[z0,z1,z2]
        a,b,c=Regression_Plane1(X,Y,Z)[0]
        fun1=lambda X: (a*X[0]+b*X[1]+c)
        #return z0
        return fun1([pressure,temperature])
            
            
            
         
    def find_pressure_index(self,pressure):#Bar
        '''
        Find index (lower bound)
        '''
        pp = 300000/self.Model.index_P
        p= int((pressure)/pp)
        if p>=self.Model.index_P:
            return self.Model.index_P-1
        else:
            return p
        
        
    def find_temperature_index(self,temperature):#K
        '''
        Find index (lower bound)
        '''
        tt = 2000/self.Model.index_T
        t=int((temperature-773)/tt)

        if t >=self.Model.index_T:
            return self.Model.index_T-1
        else:
            return t  
    
    def return_T(self):
        return self.Temperature
        

        
    
    
    def plot_calculate_velosity(self,Check_water=True,OLwater=None,WAwater=None,RIwater=None,choice=3,Phase=None,changeKprime=None):
        self.average=choice
        coee=1e5
        if Phase is not None:
            self.Phase_diagram = np.array(Phase)
            coee=1e9
            self.Pressure=self.Phase_diagram[:,0]
            self.Temperature=self.Phase_diagram[:,2]
            self.Depth = self.Phase_diagram[:,1]
        number = len(self.Phase_diagram)
        #self.progress.deleteLater()
        
        self.Vp = np.zeros(number);self.Vs = np.zeros(len(self.Vp));self.K=np.zeros(len(self.Vp));self.G=np.zeros(len(self.Vp));self.Rho=np.zeros(len(self.Vp))

        if Check_water == False:
            mineral_list=[an,ab,sp,hc,en,fs,mgts,odi,hpcen,hpcfs,di,he,cen,cats,jd,py,al,gr,mgmj,jdmj, \
                          capv,fo,fa,mgwa,fewa,mgri,feri,mgil,feil,co,mgpv,fepv,alpv,mppv,fppv,appv,mgcf,fecf, \
                          nacf,pe,wu,qtz,coes,st,an,ky,neph]            
            for number,phase in enumerate(self.Phase_diagram):

                #phase[53]=0;phase[54]=0;
                K=[];G=[];Rho=[];V=[]     
                for num,i in enumerate(mineral_list,start=3):
                    if phase[num] != 0:    
                        r,k,g,v,rho = i.EOS(self.Pressure[number]*coee,self.Temperature[number])
                        if Phase is not None:
                            V.append(phase[num]*v)
                        else:
                            V.append(phase[num])
                        K.append(k);G.append(g);Rho.append(rho);
#==============================================================================
#                         a,b,_ = phase[num]*v*np.array(i.Vp_Vs(self.Pressure[number]*coee,self.Temperature[number]))
#                         phase[53]+=a
#                         phase[54]+=b
#==============================================================================
                        
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
                if self.average == 7:
                    self.Vp[number],self.Vs[number],self.Rho[number] = phase[53],phase[54],KGRho.ReturnRho()
                
            #print ('calcualted')
        else:
            OL=OL_(olivine)
            WA=WA_(wadsleyte)
            RI=RI_(ringwoodite) 
            if changeKprime is not None:
                OL.change_prime(changeKprime[0])
                WA.change_prime(changeKprime[1])
                RI.change_prime(changeKprime[2])
                
            mineral_list=[an,ab,sp,hc,en,fs,mgts,odi,hpcen,hpcfs,di,he,cen,cats,jd,py,al,gr,mgmj,jdmj, \
                          capv,fo,fa,mgwa,fewa,mgri,feri,mgil,feil,co,mgpv,fepv,alpv,mppv,fppv,appv,mgcf,fecf, \
                          nacf,pe,wu,qtz,coes,st,an,ky,neph,OL,WA,RI]        
            for number,phase in enumerate(self.Phase_diagram):    
                OL.Set_Water_Iron_Condition(OLwater[number],1-phase[24]/(phase[24]+phase[25]+1e-11))
                WA.Set_Water_Iron_Condition(WAwater[number],1-phase[26]/(phase[26]+phase[27]+1e-11))
                RI.Set_Water_Iron_Condition(RIwater[number],1-phase[28]/(phase[28]+phase[29]+1e-11))
                phase[24]=0;phase[25]=0
                phase[26]=0;phase[27]=0
                phase[28]=0;phase[29]=0
                K=[];G=[];Rho=[];V=[]
                for num,i in enumerate(mineral_list,start=3):
                    if phase[num] != 0:    
                        r,k,g,v,rho = i.EOS(self.Pressure[number]*coee,self.Temperature[number])
                        K.append(k);G.append(g);Rho.append(rho);V.append(phase[num])                       
                        a,b,c=i.Vp_Vs(self.Pressure[number]*coee,self.Temperature[number])                        
                        i.Store_Vp_Vs(a,b,phase[1],phase[num])
                KGRho = Velocity_calculator(K,G,Rho,V)
                if self.average == 1:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()   
                    self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.

                if self.average == 2:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()   
                    self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.

                if self.average == 3:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()   
                    self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.

                if self.average == 4:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()   
                    self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.

                if self.average == 5:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()   
                    self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.
                 
                if self.average == 6:
                    self.K[number],self.G[number],self.Rho[number] =  KGRho.Hashinand_Shtrikman_lower()   
                    self.Vp[number],self.Vs[number]=np.sqrt((self.K[number]+4.*self.G[number]/3.)/self.Rho[number])/1000. , np.sqrt(self.G[number]/self.Rho[number])/1000.
                
                if self.average == 7:
                    self.Vp[number],self.Vs[number],self.Rho[number] = phase[53],phase[54],KGRho.ReturnRho()
            #print ('calcualted')


            
            
    def plot_function(self):
        '''
        This plot the phase diagram over whole upper mantle.
        '''
        fig = plt.figure()
        self.ax = fig.add_subplot(111) 
        
        global_Mineral_V=dict()
        mineral_list=['P','D','T','an','ab','sp','hc','en','fs','mgts','odi','hpcen','hpcfs','di','he','cen','cats','jd','py','al','gr','mgmj','jdmj', \
                      'capv','fo','fa','mgwa','fewa','mgri','feri','mgil','feil','co','mgpv','fepv','alpv','mppv','fppv','appv','mgcf','fecf', \
                      'nacf','pe','wu','qtz','coes','st','apbo','ky','neph','OL','WA','RI','Vp','Vs']        
        self.Phase_diagram = np.array(self.Phase_diagram)
        if len(self.Phase_diagram[0])<=52:
            return 0
        for i,name in  enumerate(mineral_list, start=0):
            #print (i)
            global_Mineral_V[name]=self.Phase_diagram[:,i]
        self.global_Mineral_V = global_Mineral_V


        global_m=len(global_Mineral_V['ab'])
        Mineral_V_Plot=dict()
        name=['sp','cpx','opx','ol','plg','ak','fp','hpcpx','gt','wa','ri','pv','ppv','cf','q','coe','st','capv']
        for i in name:   
            Mineral_V_Plot[i]=np.zeros(global_m)
        Mineral_V_Plot['sp']=global_Mineral_V['sp']+global_Mineral_V['hc']
        Mineral_V_Plot['cpx']=global_Mineral_V['di']+global_Mineral_V['he']+global_Mineral_V['cen']+global_Mineral_V['cats']+global_Mineral_V['jd']
        Mineral_V_Plot['opx']=global_Mineral_V['en']+global_Mineral_V['fs']+global_Mineral_V['mgts']+global_Mineral_V['odi']
        try:
            Mineral_V_Plot['ol']=global_Mineral_V['OL']
            Mineral_V_Plot['wa']=global_Mineral_V['WA']
            Mineral_V_Plot['ri']=global_Mineral_V['RI']
        except:
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
            self.ax.text(global_D[int(50/a*global_m)], 0.5*Mineral_V_Plot['ol'][int(50/a*global_m)], 'ol', fontsize=12)
            self.ax.text(global_D[int(170/a*global_m)], 0.5*Mineral_V_Plot['ol'][int(50/a*global_m)], 'wa', fontsize=12)
            self.ax.text(global_D[int(220/a*global_m)], 0.5*Mineral_V_Plot['ol'][int(50/a*global_m)], 'ri', fontsize=12)
            #ax.text(global_D[int(280/a*global_m)], 0.5*Mineral_V_Plot['fp'][int(280/a*global_m)], 'fp', fontsize=12)
        self.ax.text(global_D[int(160/a*global_m)], 0.5*Mineral_V_Plot['gt'][int(160/a*global_m)]+line11[int(160/a*global_m)], 'gt', fontsize=12)
        #ax.text(global_D[int(280/a*global_m)], 0.5*Mineral_V_Plot['pv'][int(280/a*global_m)]+line11[int(280/a*global_m)], 'pv', fontsize=12)
        if Mineral_V_Plot['opx'][int(10/a*global_m)] >=0.05:
            self.ax.text(global_D[int(10/a*global_m)], 0.7*Mineral_V_Plot['opx'][int(40/a*global_m)]+line11[int(40/a*global_m)], 'opx', fontsize=12)
        if Mineral_V_Plot['hpcpx'][int(130/a*global_m)] >=0.02:
            self.ax.text(global_D[int(115/a*global_m)], 0.8*Mineral_V_Plot['hpcpx'][int(130/a*global_m)]+line11[int(130/a*global_m)]+Mineral_V_Plot['gt'][int(40/a*global_m)], 'hpcpx', fontsize=12) 
        self.ax.text(global_D[int(40/a*global_m)], 0.5*Mineral_V_Plot['cpx'][int(40/a*global_m)]+line11[int(40/a*global_m)]+Mineral_V_Plot['gt'][int(40/a*global_m)]+Mineral_V_Plot['opx'][int(40/a*global_m)], 'cpx', fontsize=12) 
        if Mineral_V_Plot['st'][int(150/a*global_m)]>=0.02:
            self.ax.text(global_D[int(40/a*global_m)], 1-0.7*Mineral_V_Plot['st'][int(150/a*global_m)], 'co', fontsize=12)
            self.ax.text(global_D[int(150/a*global_m)], 1-0.7*Mineral_V_Plot['st'][int(150/a*global_m)], 'st', fontsize=12)
        if Mineral_V_Plot['capv'][int(230/a*global_m)]>=0.02:
            self.ax.text(global_D[int(230/a*global_m)], 1-0.7*Mineral_V_Plot['capv'][int(230/a*global_m)], 'capv', fontsize=12)

        self.ax.set_xlim(10,global_D[-1])
        self.ax.set_ylim(0,1)    
        self.ax.set_xlabel('Depth (km)')
        self.ax.set_ylabel('Phase proportion (volume%)')
        return fig
        
