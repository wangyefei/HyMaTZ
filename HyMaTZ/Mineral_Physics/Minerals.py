# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 16:51:45 2016

@author: Fei
Update: 10/3/2016
"""
import numpy as np
from scipy.constants import gas_constant
try:
    from tools import func,dictionarize_formula
except:
    from .tools import func,dictionarize_formula




class Mineral(object):
    """
    This class use to store mineral's data and Methods
    """
    
    def __init__(self,pressure=None,temperature=None):
        self.pressure=pressure
        self.temperature=temperature  
        self.Depth=[]
        self.Vp=[]
        self.Vs=[]
        self.Rho=[]
        self.Precentage=[]
        pass

    def Formular(self):
        a=[0,0,0,0,0,0];j=0
        aa = dictionarize_formula(self.formula)
        for i in ['Ca','Al','Na','Mg','Fe','Si']:
            try:
                a[j] = aa[i];j+=1
            except:
                j+=1
        return np.array(a)
    
    def Composition(self,molar=None):
        a=[0,0,0,0,0,0];j=0
        aa = dictionarize_formula(self.formula)
        for i in ['Na','Mg','Al','Si','Ca','Fe']:
            
            try:
                a[j] = aa[i];j+=1
            except:
                j+=1
        return np.array(a)*molar
    
    def Atom_Weight1(self,molar=1):
        a=[0,0,0,0,0,0];j=0
        aa = dictionarize_formula(self.formula)
        for i in ['Ca','Al','Na','Mg','Fe','Si']:        
            try:
                a[j] = aa[i];j+=1
            except:
                j+=1
        return np.array(a)*molar        
    
    def Atom_Weight(self):
        a=np.array(self.Formular())
       # b=np.array([56.0774, 101.9612/2,61.979/2, 40.3044,71.8444,60.0843])
        b=np.array([40.078,26.98,22.989,24.3050,55.854,28.085])
        Weight_pre=(a*b)/sum(a*b)
        return Weight_pre    
    
    def Change_name_formula(self,string1,string2):
        self.name = string1
        self.formula = string2
 
    def Clear_Vp_Vs(self):
        self.Vp=[]
        self.Vs=[]   
        self.Rho=[]
        self.Depth=[]
        self.Precentage=[]

    def Store_Vp_Vs(self,Vp,Vs,Rho,Depth,Precentage=0):
        self.Vp.append(Vp)
        self.Vs.append(Vs)
        self.Rho.append(Rho)
        self.Depth.append(Depth)
        self.Precentage.append(Precentage)
    
    def Set_Store_Vp_Vs(self,Vp,Vs,Depth,Precentage):
        self.Vp=Vp
        self.Vs=Vs        
        self.Depth=Depth
        self.Precentage=Precentage        

    def set_PT(self,pressure,temperature):
        '''
        Re set Pressure and temperature
        '''
        self.pressure=pressure
        self.temperature=temperature     
        return None
    

    


        
class NAMs(object):
    '''
    This class use to store Nominally Anhydrous Minerals Methods
    If change water or iron condition, change Minerlas properties instantly,
    '''
    
    def __init__(self,water=0,iron=0):
        self.water=water;self.iron=iron
    
    
    def Set_Water_Iron_Condition(self,water,iron,random=[0,0,0]):
        self.water = water
        self.iron = iron
        self.K_0=func(self.a[0][0],[self.water,self.iron])*1e9+random[0]*self.uncertainties['err_K_0']
        self.G_0=func(self.a[1][0],[self.water,self.iron])*1e9+random[1]*self.uncertainties['err_G_0']
        self.Rho_0=func(self.a[2][0],[self.water,self.iron])
        self.V_0=(self.molar_mass/self.Rho_0)+random[2]*self.uncertainties['err_V_0']
    
     

    def GUI_Plot_function(self,P,T):
        pass

    def GUI_clean(self,P,T):
        self.Depth=[]
        self.Vp=[]
        self.Vs=[]

class NOUSE():
    
    def __init__(self):
        pass
    
    def EOS(self):
        return 0,0,0,0,0
    
    def Vp_Vs(self):
        return 0,0,0
    
NOUSE = NOUSE()   
if __name__ == "__main__":
    pass        