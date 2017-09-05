# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 12:48:20 2016

@author: Fei

Update 10/3/2016
"""
import sys
import os
import re
from fractions import Fraction
import numpy as np
from scipy.constants import gas_constant
from scipy.optimize import nnls
try:
    from  Stix2011data import (OL,WA,RI,ab,an,sp,hc,fo,fa,mgwa,fewa,mgri,feri,en,fs,mgts,odi,di\
                          ,he,cen,cats,jd,hpcen,hpcfs,mgpv,fepv\
                          ,alpv,capv,mgil,feil,co,py ,al,gr,mgmj ,jdmj,qtz,coes,st,seif,mppv,fppv,appv,pe,wu,mgcf \
                          ,fecf ,nacf,ky,neph )
except:
    from  .Stix2011data import (OL,WA,RI,ab,an,sp,hc,fo,fa,mgwa,fewa,mgri,feri,en,fs,mgts,odi,di\
                          ,he,cen,cats,jd,hpcen,hpcfs,mgpv,fepv\
                          ,alpv,capv,mgil,feil,co,py ,al,gr,mgmj ,jdmj,qtz,coes,st,seif,mppv,fppv,appv,pe,wu,mgcf \
                          ,fecf ,nacf,ky,neph )

'''
Solid Solution contain each solid solution's minerals, Modify from Burnman please also see Burnman.
Minerals from Stixrude & Lithgow-Bertelloni 2005, 2011 and references therein,
''' 


class SolidSolution(object):
    """
    This is the base class for all solid solutions.
    
    """
    def __init__(self,molar_fractions=None):
        self.state = False        
        if molar_fractions == None:
            self.molar_fractions = np.zeros(len(self.endmembers))
            self.molar_fractions[0]=1
        else:
            self.molar_fractions=molar_fractions
        return None
        
    def GUI_Plot_function(self):            
        self.Depth = self.endmembers[0][0].Depth
        for mineral in self.endmembers:
            if len(mineral[0].Depth) >= len(self.Depth):
                self.Depth = mineral[0].Depth


        for mineral in self.endmembers:
            Vp = np.zeros(len(self.Depth))
            Vs = np.zeros(len(self.Depth))
            Precentage = np.zeros(len(self.Depth))
            if len(mineral[0].Depth) == 0:
                mineral[0].Set_Store_Vp_Vs(Vp,Vs,self.Depth,Precentage) 
            elif len(mineral[0].Depth) == len(self.Depth):
                pass                
            else:
                for num,depth in enumerate(self.Depth):                
                    if depth == mineral[0].Depth[0]:
                        a = num
                    if depth == mineral[0].Depth[-1]:
                        b = num
                Vp[a:b+1] = np.array(mineral[0].Vp) 
                Vs[a:b+1] = np.array(mineral[0].Vs) 
                Precentage[a:b+1]=np.array(mineral[0].Precentage)
                mineral[0].Set_Store_Vp_Vs(Vp,Vs,self.Depth,Precentage)

        self.Vp = np.zeros(len(self.Depth))
        self.Vs = np.zeros(len(self.Depth))
        self.Precentage = np.zeros(len(self.Depth))  
        sum1 = np.zeros(len(self.Depth))   
        for mineral in self.endmembers:
            sum1 += mineral[0].Precentage
        for mineral in self.endmembers:
            self.Vp += mineral[0].Vp * (mineral[0].Precentage/sum1)
            self.Vs += mineral[0].Vs * (mineral[0].Precentage/sum1)
            

    def GUI_clean(self):
        self.Depth=[]
        self.Vp=[]
        self.Vs=[]
            
        
    def Composition(self,molar_fractions=None):
        if molar_fractions is None:
            return [0,0,0,0,0,0,]
        else:
            a = np.array([0.5,1,0.5,1,1,1])
            molar = np.zeros(6)
            for i in range(len(self.composition)):
                molar += np.array(self.composition[i])*a*molar_fractions[i]  
            return molar
        

    
    def Volume(self,Pressure,Temperature,Composition,dp=1e-1):
        return sum(self.Gibbs_energy(Pressure+dp,Temperature,Composition)*Composition-self.Gibbs_energy(Pressure-dp,Temperature,Composition)*Composition)/2/dp        
        
        
    def change_state(self,state):
        self.state = state
        
    def Set_Molar_fractions(self,molar_fractions=None):
        self.molar_fractions=molar_fractions
        return None
        
    
    def Endmembers_Matrix(self):
        w= len(self.endmembers)
        Matrix = np.array([self.endmembers[x][0].Formular() for x in range(w)])
        return Matrix     
    


    def Solve_Molar_fraction(self,weight_fraction=[0,0,0,0,0,0]):
        if sum(weight_fraction)==0:
            return np.zeros(len(self.endmembers))
        a=np.zeros((len(self.endmembers),6))
        for i in range(len(self.endmembers)):
            a[i]=self.endmembers[i][0].Atom_Weight()
        a=a.transpose()
        x = nnls(a, np.array(weight_fraction))[0]
        for i in range(len(x)):
            check=True
            if len(x) != len(self.endmembers):continue
            for j in range(len(self.endmembers)):
                if x[j] < 0.0 :
                    check=False
        if check==True:
            return np.array([i/sum(x) for i in x])
        else:
            print ('bad end')
            return x
        
    def test(self):
        a=np.zeros((len(self.endmembers),6))
        for i in range(len(self.endmembers)):
            a[i]=self.endmembers[i][0].Atom_Weight()
        a=a.transpose()
        return a                


        



class c2c_pyroxene(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'C2/c pyroxene'
        self.type = 'ideal'
        self.endmembers = [
            [hpcen, '[Mg]2Si2O6'], [hpcfs, '[Fe]2Si2O6']]
        self.energy_interaction = [[0.0]]
        self.composition = [[0, 2.0, 0, 2.0, 0, 0],[0, 0, 0, 2.0, 0, 2.0]]
        

        SolidSolution.__init__(self, molar_fractions)


class ca_ferrite_structured_phase(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'calcium ferrite structured phase'
        self.type = 'ideal'
        self.endmembers = [[mgcf, '[Mg]Al[Al]O4'], [
                           fecf, '[Fe]Al[Al]O4'], [nacf, '[Na]Al[Si]O4']]
        self.energy_interaction = [[0.0, 0.0], [0.0]]
        self.composition = [[0, 1.0, 2.0, 0, 0, 0],[0, 0, 2.0, 0, 0, 1.0],[1.0, 0, 2.0, 0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class clinopyroxene(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'clinopyroxene'
        self.type = 'symmetric'
        self.endmembers = [[di, '[Ca][Mg][Si]2O6'], [he, '[Ca][Fe][Si]2O6'], [
                           cen, '[Mg][Mg][Si]2O6'], [cats, '[Ca][Al][Si1/2Al1/2]2O6'], [jd, '[Na][Al][Si]2O6']]
        self.energy_interaction = [
            [0., 24.74e3, 26.e3, 24.3e3], [24.74e3, 0., 0.e3], [60.53136e3, 0.0], [10.e3]]
        self.composition = [[0, 1.0, 0, 2.0, 1.0, 0],[0, 0, 0, 2.0, 1.0, 1.0],[0, 2.0, 0, 2.0, 0, 0],[0, 0, 2.0, 1.0, 1.0, 0],[1.0, 0, 1.0, 2.0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class garnet(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'garnet'
        self.type = 'symmetric'
        self.endmembers = [[py, '[Mg]3[Al][Al]Si3O12'], [al, '[Fe]3[Al][Al]Si3O12'], [
                                   gr, '[Ca]3[Al][Al]Si3O12'], [mgmj, '[Mg]3[Mg][Si]Si3O12'], [jdmj, '[Na2/3Al1/3]3[Al][Si]Si3O12']]
        self.energy_interaction = [
            [0.0, 30.e3, 21.20278e3, 0.0], [0.0, 0.0, 0.0], [57.77596e3, 0.0], [0.0]]
        self.composition = [[0, 3.0, 2.0, 3.0, 0, 0],[0, 0, 2.0, 3.0, 0, 3.0],[0, 0, 2.0, 3.0, 3.0, 0],[0, 4.0, 0, 4.0, 0, 0],[2.0, 0, 2.0, 4.0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class akimotoite(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'akimotoite/ilmenite'
        self.type = 'symmetric'
        self.endmembers = [[mgil, '[Mg][Si]O3'], [
                           feil, '[Fe][Si]O3'], [co, '[Al][Al]O3']]
        self.energy_interaction = [[0.0, 66.e3], [66.e3]]
        self.composition = [[0, 1.0, 0, 1.0, 0, 0],[0, 0, 0, 1.0, 0, 1.0],[0, 0, 2.0, 0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class ferropericlase(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'magnesiowustite/ferropericlase'
        self.type = 'symmetric'
        self.endmembers = [[pe, '[Mg]O'], [wu, '[Fe]O']]
        self.energy_interaction = [[13.e3]]
        self.composition = [[0, 1.0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 1.0]]
        SolidSolution.__init__(self, molar_fractions)





class orthopyroxene(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'orthopyroxene'
        self.type = 'symmetric'
        self.endmembers = [[en, '[Mg][Mg][Si]SiO6'], [fs, '[Fe][Fe][Si]SiO6'], [
                           mgts, '[Mg][Al][Al]SiO6'], [odi, '[Ca][Mg][Si]SiO6']]
        self.energy_interaction = [
            [0.0, 0.0, 32.11352e3], [0.0, 0.0], [48.35316e3]]
        self.composition = [[0, 2.0, 0, 2.0, 0, 0],[0, 0, 0, 2.0, 0, 2.0],[0, 1.0, 2.0, 1.0, 0, 0],[0, 1.0, 0, 2.0, 1.0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class plagioclase(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'plagioclase'
        self.type = 'symmetric'
        self.endmembers = [
            [an, '[Ca][Al]2Si2O8'], [ab, '[Na][Al1/2Si1/2]2Si2O8']]
        self.energy_interaction = [[26.0e3]]
        self.composition = [[0, 0, 2.0, 2.0, 1.0, 0],[1.0, 0, 1.0, 3.0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class post_perovskite(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'post-perovskite/bridgmanite'
        self.type = 'symmetric'
        self.endmembers = [[mppv, '[Mg][Si]O3'], [
                           fppv, '[Fe][Si]O3'], [appv, '[Al][Al]O3']]
        self.energy_interaction = [[0.0, 60.0e3], [0.0]]
        self.composition = [[0, 1.0, 0, 1.0, 0, 0],[0, 0, 0, 1.0, 0, 1.0],[0, 0, 2.0, 0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)


class mg_fe_perovskite(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'magnesium silicate perovskite/bridgmanite'
        self.type = 'symmetric'
        self.endmembers = [[mgpv, '[Mg][Si]O3'], [
                           fepv, '[Fe][Si]O3'], [alpv, '[Al][Al]O3']]
        self.energy_interaction = [[0.0, 116.0e3], [0.0]]
        self.composition = [[0, 1.0, 0, 1.0, 0, 0],[0, 0, 0, 1.0, 0, 1.0],[0, 0, 2.0, 0, 0, 0]]
        SolidSolution.__init__(self, molar_fractions)
        
        
class mg_fe_aluminous_spinel(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'spinel-hercynite binary, fixed order'
        self.type = 'symmetric'
        self.endmembers = [[sp, '[Mg3/4Al1/4]4[Al7/8Mg1/8]8O16'], [
                                   hc, '[Fe3/4Al1/4]4[Al7/8Fe1/8]8O16']]
        self.energy_interaction = [[5.87646e3]]
        self.composition = [[0, 3.0, 7.0, 0, 0, 0],[0, 0, 7.0, 0, 0, 3.0]]
        SolidSolution.__init__(self, molar_fractions)


    

class mg_fe_olivine(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'olivine'
        self.type = 'symmetric'
        self.endmembers = [[
            fo, '[Mg]2SiO4'], [fa, '[Fe]2SiO4']]
        self.energy_interaction = [[7.81322e3]]
        self.composition = [[0, 2.0, 0, 1.0, 0, 0],[0, 0, 0, 1.0, 0, 2.0]]
        SolidSolution.__init__(self, molar_fractions)

    def Gibbs_water(self,P,T,mole_fractions,water_WT=0):
#==============================================================================
#         XOH  = water_WT/6.6
#         a = (1-XOH)**4
#         b = (1-0.5*XOH)**4 
#         return gas_constant*T*np.log(a/b)/1000.        
#         pass
#==============================================================================
        water_molr_fraction = water_WT/3.3 *0.5
        XOH = water_molr_fraction/2.
        return gas_constant*T*np.log((1-XOH)**0.5)/1000.
    
    def Gibbs(self,P,T,mole_fractions,Tr=300,water_WT=0):
        a=self.Gibbs_energy(P,T,mole_fractions)
        if sum(mole_fractions)==0:
            return 0
        else:
            mole_fractions/=sum(np.array(mole_fractions))
            return sum(a*mole_fractions)+self.Gibbs_water(P,T,mole_fractions,water_WT)  

    
class mg_fe_wadsleyite(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'wadsleyite'
        self.type = 'symmetric'
        self.endmembers = [
            [mgwa, '[Mg]2SiO4'], [fewa, '[Fe]2SiO4']]
        self.energy_interaction = [[16.74718e3]]
        self.composition = [[0, 2.0, 0, 1.0, 0, 0],[0, 0, 0, 1.0, 0, 2.0]]
        SolidSolution.__init__(self, molar_fractions)
        
    def Gibbs_water(self,P,T,mole_fractions,water_WT=0):
        water_molr_fraction = water_WT/3.3 *0.5
        XOH = water_molr_fraction/2.
        return gas_constant*T*np.log((1-XOH)**0.5)/1000.
    
    def Gibbs(self,P,T,mole_fractions,Tr=300,water_WT=0):
        a=self.Gibbs_energy(P,T,mole_fractions)
        if sum(mole_fractions)==0:
            return 0
        else:
            mole_fractions/=sum(np.array(mole_fractions))
            return sum(a*mole_fractions)+self.Gibbs_water(P,T,mole_fractions,water_WT)        

        
class mg_fe_ringwoodite(SolidSolution):#

    def __init__(self, molar_fractions=None):
        self.name = 'ringwoodite'
        self.type = 'symmetric'
        self.endmembers = [
            [mgri, '[Mg]2SiO4'], [feri, '[Fe]2SiO4']]
        self.energy_interaction = [[9.34084e3]]
        self.composition = [[0, 2.0, 0, 1.0, 0, 0],[0, 0, 0, 1.0, 0, 2.0]]
        SolidSolution.__init__(self, molar_fractions)

    def Gibbs_water(self,P,T,mole_fractions,water_WT=0):
        water_molr_fraction = water_WT/3.3 *0.5
        XOH = water_molr_fraction/2.
        return gas_constant*T*np.log((1-XOH)**0.5)/1000.
    
    def Gibbs(self,P,T,mole_fractions,Tr=300,water_WT=0):
        a=self.Gibbs_energy(P,T,mole_fractions)
        if sum(mole_fractions)==0:
            return 0
        else:
            mole_fractions/=sum(np.array(mole_fractions))
            return sum(a*mole_fractions)+self.Gibbs_water(P,T,mole_fractions,water_WT)           



class OL_water(SolidSolution):
    
    def __init__(self, OL):
        self.name = 'olivine'
        self.endmembers = [[OL,'(Mg,Fe)2SiO4']]

class WA_water(SolidSolution):
    
    def __init__(self, WA):
        self.name = 'wadsleyite'
        self.endmembers = [[WA,'(Mg,Fe)2SiO4']]
        
class RI_water(SolidSolution):
    
    def __init__(self, RI):
        self.name = 'ringwoodite'
        self.endmembers = [[RI,'(Mg,Fe)2SiO4'] ]   



# Solid solution aliases
c2c = c2c_pyroxene()   #2  C2/c
CF = ca_ferrite_structured_phase()#2
Cpx = clinopyroxene() #5
Gt = garnet()#5
Gt_maj = garnet()#5
Aki = akimotoite()#3
Wus = ferropericlase()#2
O = mg_fe_olivine()
Opx = orthopyroxene()#4
Pl = plagioclase()#2
ppv = post_perovskite()#3
Ppv = post_perovskite()
Pv = mg_fe_perovskite()#3
Ring = mg_fe_ringwoodite()
Sp = mg_fe_aluminous_spinel()
Wad = mg_fe_wadsleyite()
OLwater = OL_water(OL)
WAwater = WA_water(WA)
RIwater = RI_water(RI)

if __name__ == "__main__":
    pass

