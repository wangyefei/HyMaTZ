# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 20:33:49 2017
Methods is the same as reference and the colde in EOS, but stored data based on different EOS
@author: Fei
"""

from __future__ import print_function
import numpy as np
from scipy import integrate
from scipy.optimize import newton,bisect
from scipy.constants import gas_constant


try:
    from tools import dictionarize_formula,formula_mass,atomic_masses
    from Minerals import Mineral,NAMs
    from regression import olivine,wadsleyte,ringwoodite
    
except:
    from .tools import (dictionarize_formula,formula_mass,atomic_masses)
    from .Minerals import Mineral,NAMs
    from .regression import olivine,wadsleyte,ringwoodite
    


class EOS(Mineral,NAMs):
    '''
    Basic class 
    '''
    def __init__(self):
        super(EOS,self).__init__()
        self.name = 'input name'
        self.formula = 'input formula'
        
        self.n = sum(dictionarize_formula(self.formula).values())
        self.molar_mass = formula_mass(dictionarize_formula(self.formula),atomic_masses)
        
        self.HPparams = {
            #'H_0': 0,
            'S_0(J/K)': 0,
            'V_0(cm³/mol)': 0,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 0,
            'K(Pa)': 0,
            'G(Pa)':0,
            "k'": 0,
            "G'":0,
            "K''": 0,
            }
        
        self.SLBparams = {
            #'F_0': 0,
            'V_0(cm³/mol)': 0,
            'K(Pa)': 0,
            "k'": 0,
            'Debye_0(K)': 0,
            'grueneisen_0': 0,
            'q_0': 0,
            'G(Pa)': 0,
            "G'": 0,
            'eta_s_0': 0,
            }        

        self.BM2rdparams = {
            #'H_0': 0,
            #'S_0(J/K)':0,
            'V_0(cm³/mol)': 0,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 0,
            'K(Pa)': 0,
            'G(Pa)':0,
            "k'": 0,
            "G'":0,
            "K''": 0,
            'dK/dT':0,
            'dG/dT':0,
            }

        self.BM3rdparams = {
            #'H_0': 0,
            #'S_0(J/K)': 0,
            'V_0(cm³/mol)': 0,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 0,
            'K(Pa)': 0,
            'G(Pa)':0,
            "k'": 0,
            "G'":0,
            "K''": 0,
            'dK/dT':0,
            'dG/dT':0,
            }

        self.BMparams = {
            #'H_0': 0,
            #'S_0(J/K)': 0,
            'V_0(cm³/mol)': 0,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 0,
            'K(Pa)': 0,
            'G(Pa)':0,
            "k'": 0,
            "G'":0,
            "K''": 0,
            'dK/dT':0,
            'dG/dT':0,
            } 
        
    def Print_params(self):
        line1='name;'+'formula;'+'n;'+'molar mass;'
        line2=self.name+';'+self.formula+';'+str(self.n)+';'+str(self.molar_mass)+';'
        
        line3='';line4=''
        for key,value in self.SLBparams.items():
            line3+=str(key)+';'
            line4+=str(value)+';'
        
        line5='';line6=''
        for key,value in self.HPparams.items():
            line5+=str(key)+';'
            line6+=str(value)+';'   
        
        line7='';line8=''
        for key,value in self.BM2rdparams.items():
            line7+=str(key)+';'
            line8+=str(value)+';'
        
        line9='';line10=''
        for key,value in self.BM3rdparams.items():
            line9+=str(key)+';'
            line10+=str(value)+';'
        
        line11='';line12=''
        for key,value in self.BMparams.items():
            line11+=str(key)+';'
            line12+=str(value)+';'
            
        return line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12
    
    
    def Restart(self,lines):
        line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12 = lines
        line2 = line2.split(';')
        self.name = line2[0];self.formula=line2[1];self.n=float(line2[2]);self.molar_mass=float(line2[3])
        
        line3 = line3.split(';')[:-1]; line4 = line4.split(';')[:-1]
        print (line4,line3)
        for num,key in enumerate(line3):
            self.SLBparams[key]=float(line4[num])

        line5 = line5.split(';')[:-1]; line6 = line6.split(';')[:-1]
        for num,key in enumerate(line5):
            self.HPparams[key]=float(line6[num])      

        line7 = line7.split(';')[:-1]; line8 = line8.split(';')[:-1]
        for num,key in enumerate(line7):
            self.BM2rdparams[key]=float(line8[num]  ) 
            
        line9 = line9.split(';')[:-1]; line10 = line10.split(';')[:-1]
        for num,key in enumerate(line9):
            self.BM3rdparams[key]=float(line10[num]    )           

        line11 = line11.split(';')[:-1]; line12 = line12.split(';')[:-1]
        for num,key in enumerate(line11):
            self.BMparams[key]=float(line12[num] )
        print ('Restart Done')
        return None
    
    def Change_SLB(self,list1):
        for num,key in enumerate(self.SLBparams.keys()):
            self.SLBparams[key] = list1[num]
            
    def Change_HP(self,list1):
        for num,key in enumerate(self.HPparams.keys()):
            self.HPparams[key] = list1[num]

    def Change_BM2rd(self,list1):
        for num,key in enumerate(self.BM2rdparams.keys()):
            self.BM2rdparams[key] = list1[num]

    def Change_BM3rd(self,list1):
        for num,key in enumerate(self.BM3rdparams.keys()):
            self.BM3rdparams[key] = list1[num]

    def Change_BM(self,list1):
        for num,key in enumerate(self.BMparams.keys()):
            self.BMparams[key] = list1[num]            
            
    def Change_Name_formula(self,name,formula):
        self.name=name  
        self.formula = formula
        self.n = sum(dictionarize_formula(self.formula).values())
        self.molar_mass=formula_mass(dictionarize_formula(self.formula),atomic_masses)  
        
    def Change_Name(self,name):
        self.name=name   
    
    def Change_formula(self,formula):
        self.formula = formula
        self.n = sum(dictionarize_formula(self.formula).values())
        self.molar_mass=formula_mass(dictionarize_formula(self.formula),atomic_masses)  
        
    '''
    TEOS from Holand and Powell (2011)
    '''
    def HP_K_abs1(self):#K0，K1,K2 -> abc
        a=(1+self.HPparams["k'"])/(1+self.HPparams["k'"]+self.HPparams['K(Pa)']*self.HPparams["K''"])
        b=(self.HPparams["k'"]/self.HPparams['K(Pa)'])-(self.HPparams["K''"]/(1+self.HPparams["k'"]))
        c=(1+self.HPparams["k'"]+self.HPparams['K(Pa)']*self.HPparams["K''"])/(self.HPparams["k'"]**2+self.HPparams["k'"]-self.HPparams['K(Pa)']*(self.HPparams["K''"]**2))
        return a,b,c       

    def HP_VP0(self,Pressure):#return VP0
        a,b,c=self.HP_K_abs1()
        VP0=(1-a*(1-(1+b*Pressure)**(-c)))*self.HPparams['V_0(cm³/mol)']
        return VP0    

    def HP_Pth(self,Temperature,Tr=300):
        a,b,c=self.HP_K_abs1()
        TD = 10636/(self.HPparams['S_0(J/K)']/self.n+6.44)
        u=TD/Temperature;u0=TD/Tr
        segma=(u0*u0*np.exp(u0))/(np.exp(u0)-1)**2
        Pth=(float(self.HPparams['a_0(K-1)'])*self.HPparams['K(Pa)']*TD/segma)*(1/(np.exp(u)-1)-1/(np.exp(u0)-1))
        return Pth #return Pth    
    
    def HP_VPT(self,Pressure,Temperature,Tr=298.):
        a,b,c=self.HP_K_abs1()    
        Pth=self.HP_Pth(Temperature,Tr) 
        VPT=(1-a*(1-(1+b*(Pressure-Pth))**(-c)))*self.HPparams['V_0(cm³/mol)']
        return VPT 

    def HP_K(self,Pressure,Temperature,Tr=298.):
        a,b,c=self.HP_K_abs1()    
        Pth=self.HP_Pth(Temperature,Tr) 
        K=self.HPparams['K(Pa)']*(1+b*(Pressure-Pth))*(a+(1-a)*(1+b*(Pressure-Pth))**c)  
        return K

    def HP_G(self,Pressure,Temperature,Tr=298.,return_Rho=False):
        '''
        Shear modulus from Connolly & Kerrick (2002) and Connolly (2005)
        '''

        V=self.HP_VPT(Pressure,Temperature,Tr) 
        r=self.HPparams['V_0(cm³/mol)']/V               
        f = 0.5 * (pow(r, 2. / 3.) - 1.0)               
        #G=self.params['G(Pa)'] + (Temperature-Tr)*self.params['Gprime_0T'] + Pressure*self.params["G'"]
        G=(1+2.*f)**(5./2.)*(self.HPparams['G(Pa)']+(3.*self.HPparams['K(Pa)']*self.HPparams["G'"]-5.*self.HPparams['G(Pa)'])*f+(6.*self.HPparams['K(Pa)']*self.HPparams["G'"]-24.*self.HPparams['K(Pa)']-14.*self.HPparams['G(Pa)']+4.5*self.HPparams['K(Pa)']*self.HPparams["k'"])*f*f)#-ets*(EthVT)*rho        
        Rho = Rho=self.molar_mass/V
        if return_Rho:
            return G,Rho
        else:
            return G
        
    
    def HP_Vp_Vs(self,Pressure,Temperature,Tr=298):
        '''
        Returns unit:[km/s]`
        Returns unit:[km/s]`
        Returns unit:[kg/m3]`
        '''
        K=self.HP_K(Pressure,Temperature)
        G,Rho = self.HP_G(Pressure,Temperature,Tr,return_Rho=True)
        Vp=np.sqrt((K+4.*G/3.)/Rho)/1000.
        Vs=np.sqrt(G/Rho)/1000.;#print Vs
        return Vp,Vs,Rho/1000.,K,G 
        
    '''
    Mie–Gr¨uneisen–Debye formalism from Stixrude & Lithgow-Bertelloni (2005)
    '''
    def Stix_Compressive(self,P,T,Tr=300.,term='isotropic'):
        '''
        Given pressure (P) and temperature (T),return parameters
        Units is based on input, default is SI units
        r = V0/VPT at PT      No units   
        Units can changed based on input
        '''        
        R=gas_constant
        #eq 25，26 2011
        Debye_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)  # from wiki
        aii=6.*self.SLBparams['grueneisen_0']; #Eq47
        aiikk=-12.*self.SLBparams['grueneisen_0']+36.*(self.SLBparams['grueneisen_0']**2)-18.*self.SLBparams['q_0']*self.SLBparams['grueneisen_0']
        
        def fun(r): #r=rho/rho0
            #K2=(-1./self.params['K(Pa)'])*((3-self.params["k'"])*(4-self.params["k'"])+35./9.); # Angel 2000
            f=0.5*((r)**(2./3.)-1) #eq 20 2011
            Gru=(1./6.)/(1+aii*f+0.5*aiikk*f*f)*(2*f+1)*(aii+aiikk*f)  #  eq41, 44 2005               
            Debye=abs(np.sqrt(1.+aii*f+0.5*aiikk*f*f)*self.SLBparams['Debye_0(K)']) #eq24 2011
                       
            #use scipy to integrae Debye function 
            EthVT =(9.*self.n*R*T)  * ((Debye/T)**(-3))*integrate.quad(Debye_int_func,0.,Debye/T)[0]
            EthVT0=(9.*self.n*R*Tr) * ((Debye/Tr)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]
                        
            Pth=(EthVT-EthVT0)*(Gru*r/self.SLBparams['V_0(cm³/mol)'])  #
            
            if self.SLBparams['V_0(cm³/mol)']!=0.0:
                f=0.5*(r**(2./3.)-1.)
                a=9.*self.SLBparams['K(Pa)'];b=27.*self.SLBparams['K(Pa)']*(self.SLBparams["k'"]-4.);c=0#81.0*(self.params['K(Pa)']*K2+self.params["k'"]*(self.params["k'"]-7.)+143./9.); # eq 21,22 2011
                P0=((1./3.)*((1+2*f)**(5./2.)))*(a*f+0.5*b*f*f+1./6.*c*f*f*f)

                return P-P0-Pth#-1e5
            else:
                return 10.0
             
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)
            
        return r,aii,aiikk
    

    def Stix_thermal_pressure(self,Pressure,Temperature,Tr=300):
        r,aii,aiikk= self.Stix_Compressive(Pressure,Temperature,Tr)
        f=0.5*(r**(2./3.)-1);
        Gru=(1./6.)/(1+aii*f+0.5*aiikk*f*f)*(2*f+1)*(aii+aiikk*f)
        
        Debye=abs(np.sqrt(1+aii*f+0.5*aiikk*f*f))*self.SLBparams['Debye_0(K)']
        Debye_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)
        
        EthVT =(9.*self.n*gas_constant*Temperature)  * ((Debye/Temperature)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Temperature)[0]
        EthVT0=(9.*self.n*gas_constant*Tr) * ((Debye/Tr)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]
        Pth=(EthVT-EthVT0)*(Gru*r/self.SLBparams['V_0(cm³/mol)'])
        return Pth
        
        
        
        
    
    def Stix_EOS(self,P,T,Tr=300.5,term='isotropic',return_control = 'needed'):
        '''
        Given pressure (P) and temperature (T),return parameters
        Units is based on input, default is SI units
        r = V0/VPT at PT      No units
        K = bulk modulus at PT     Pa
        G = Shear modulus at PT    Pa
        V = Volume at PT
        Rho = Density at PT   
        a0 = Thermal expansion 
        Debye = Debye temperature K
        Units is based on input, default is 
        '''
        if P<=0:
            P==1
        R=gas_constant
        #eq 25，26 2011
        Debye_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)  # from wiki       
        r,aii,aiikk= self.Stix_Compressive(P,T,Tr)
        
        f=0.5*(r**(2./3.)-1.);
        rho=r/(self.SLBparams['V_0(cm³/mol)'])
        V=self.SLBparams['V_0(cm³/mol)']/r
        Rho=self.molar_mass/V
        Debye=abs(np.sqrt(1+aii*f+0.5*aiikk*f*f))*self.SLBparams['Debye_0(K)']
        
        EthVT=(9.*self.n*R*T)  * ((Debye/T)**(-3))*integrate.quad(Debye_int_func,0.,Debye/T)[0]
        EthVT -=(9.*self.n*R*Tr)  * ((Debye/Tr)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]
        
        CVT = 3.0 * self.n* R * (4.0 * 3*integrate.quad(Debye_int_func,0.,Debye/T)[0]/(Debye/T)**(3.) - 3.0 * (Debye/T) / (np.exp(Debye/T) - 1.0))
        CVTr = 3.0 * self.n * R * (4.0 * 3*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]/(Debye/Tr)**(3.) - 3.0 * (Debye/Tr) / (np.exp(Debye/Tr) - 1.0))
        
        Gru=(1./6.)/(1+aii*f+0.5*aiikk*f*f)*(2*f+1)*(aii+aiikk*f) 

        a2s=-2.*self.SLBparams['grueneisen_0'] - 2.*self.SLBparams['eta_s_0']  #EQ 47
        ets=-Gru-0.5/(1.+aii*f+0.5*aiikk*f*f)*((2.*f+1.)**2.)*a2s #EQ 46
        
        q=1./9.*(18.*Gru-6.-1./2./ (1.+aii*f+0.5*aiikk*f*f) * (2.*f+1.)*(2.*f+1.)*aiikk/Gru)
        K=(1+2.*f)**(5./2.)*(self.SLBparams['K(Pa)']+(3.*self.SLBparams['K(Pa)']*self.SLBparams["k'"]-5.*self.SLBparams['K(Pa)'])*f+(27./2.)*(self.SLBparams['K(Pa)']*self.SLBparams["k'"]-4.*self.SLBparams['K(Pa)'])*f*f)+(Gru+1-q)*Gru*rho*(EthVT)-Gru*Gru*rho*(CVT*T-CVTr*Tr)
        a0= Gru* CVT/(K*V)  
        G=(1+2.*f)**(5./2.)*(self.SLBparams['K(Pa)']+(3.*self.SLBparams['K(Pa)']*self.SLBparams["G'"]-5.*self.SLBparams['K(Pa)'])*f+(6.*self.SLBparams['K(Pa)']*self.SLBparams["G'"]-24.*self.SLBparams['K(Pa)']-14.*self.SLBparams['K(Pa)']+4.5*self.SLBparams['K(Pa)']*self.SLBparams["k'"])*f*f)-ets*(EthVT)*rho        
        if return_control == 'needed':
            return r,K,G,V,Rho
        else:
            return r,K,G,V,Rho,CVT,Gru,a0,Debye
        
    def Stix_Vp_Vs(self,pressure,temperature):
        '''
        Returns unit:[km/s]`
        Returns unit:[km/s]`
        Returns unit:[kg/m3]`
        '''
        r,K,G,V,Rho=self.Stix_EOS(pressure,temperature)
        Vp=np.sqrt((K+4.*G/3.)/Rho)/1000.
        Vs=np.sqrt(G/Rho)/1000.;#print Vs
        return Vp,Vs,Rho/1000.,K,G         

    
    '''
    B_M function  Ahmad & Alkammash (2012)
    '''
    def B_M_2rd(self,pressure):
        def fun(r):    
            return pressure-1.5*self.BM2rdparams['K(Pa)']*(r**(7./3.)-r**(5./3.))
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)    
        return r
    
    def B_M_3rd(self,pressure):
        def fun(r):    
            return pressure-1.5*self.BM3rdparams['K(Pa)']*(r**(7./3.)-r**(5./3.))*(1.+0.75*(r**(2./3.)-1.)*(self.BM3rdparams["k'"]-4.))
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)          
        return r
    
    def B_M_4rd(self,pressure):
        B1 = (self.BMparams['K(Pa)']*self.BMparams["K''"])    + (self.BMparams["k'"]-4.)*(self.BMparams["k'"]-5.)     + (59./9.)
        B2 = (3.*self.BMparams['K(Pa)']*self.BMparams["K''"]) + (self.BMparams["k'"]-4.)*(3.*self.BMparams["k'"]-13.) + (129./9.)
        B3 = (3.*self.BMparams['K(Pa)']*self.BMparams["K''"]) + (self.BMparams["k'"]-4.)*(3.*self.BMparams["k'"]-11.) + (105./9.)
        B4 = (self.BMparams['K(Pa)']*self.BMparams["K''"])    + (self.BMparams["k'"]-4.)*(self.BMparams["k'"]-3.)     + (35./9.)
        
        def fun(r):
            return pressure- (9.*self.BMparams['K(Pa)']/16.)*(-B1*r**(-5./3.)+B2*r**(-7./3.)-B3*r**(-3.)+B4*r**(-11./3.))
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)          
        return r,B1,B2,B3,B4  
    
   
    '''
    Thermal Pressure
    '''
    def Thermal_expansion_pressure(self,Temperature,Tr=300,order=3):
        if order ==2:
            if isinstance(self.BM2rdparams['a_0(K-1)'], float):
                Pth = (Temperature-Tr)*self.BM2rdparams['a_0(K-1)']*self.BM2rdparams['K(Pa)']
            if isinstance(self.BM2rdparams['a_0(K-1)'], str):
                T=Temperature
                a_T = eval(self.BM2rdparams['a_0(K-1)'])
                T=Tr
                a_Tr = eval(self.BM2rdparams['a_0(K-1)'])
                Pth=Temperature*a_T*self.BM2rdparams['K(Pa)']-Tr*a_Tr*self.BM2rdparams['K(Pa)']
        if order ==3:
            if isinstance(self.BM3rdparams['a_0(K-1)'], float):
                Pth = (Temperature-Tr)*self.BM3rdparams['a_0(K-1)']*self.BM3rdparams['K(Pa)']
            if isinstance(self.BM3rdparams['a_0(K-1)'], str):
                T=Temperature
                a_T = eval(self.BM3rdparams['a_0(K-1)'])
                T=Tr
                a_Tr = eval(self.BM3rdparams['a_0(K-1)'])
                Pth=Temperature*a_T*self.BM3rdparams['K(Pa)']-Tr*a_Tr*self.BM3rdparams['K(Pa)']
        if order ==4:
            if isinstance(self.BMparams['a_0(K-1)'], float):
                Pth = (Temperature-Tr)*self.BMparams['a_0(K-1)']*self.BMparams['K(Pa)']
            if isinstance(self.BMparams['a_0(K-1)'], str):
                T=Temperature
                a_T = eval(self.BMparams['a_0(K-1)'])
                T=Tr
                a_Tr = eval(self.BMparams['a_0(K-1)'])
                Pth=Temperature*a_T*self.BMparams['K(Pa)']-Tr*a_Tr*self.BMparams['K(Pa)']        
        return Pth
    
  
    def Thermal_Vp_Vs(self,Pressure,Temperature,Tr=300,order=3):
        Pth = self.Thermal_expansion_pressure(Temperature,Tr,order)
        K_0 = self.BMparams['K(Pa)']-((Temperature-Tr)*self.BMparams['dK/dT'] )
        G_0 = self.BMparams['G(Pa)']-((Temperature-Tr)*self.BMparams['dG/dT'] )

        if order ==2:
            K_0 = self.BM2rdparams['K(Pa)']-((Temperature-Tr)*self.BM2rdparams['dK/dT'] )
            G_0 = self.BM2rdparams['G(Pa)']-((Temperature-Tr)*self.BM2rdparams['dG/dT'] )
            r = self.B_M_2rd(Pressure-Pth);V=self.BMparams['V_0(cm³/mol)']/r 
            f = 0.5 * (pow(r, 2. / 3.) - 1.0)
            r = 1./r
            K = 0.5*K_0 *(7.*pow(r, -7./3.)-5.*pow(r, -5. / 3.))  
        
        elif order==3:
            K_0 = self.BM3rdparams['K(Pa)']-((Temperature-Tr)*self.BM3rdparams['dK/dT'] )
            G_0 = self.BM3rdparams['G(Pa)']-((Temperature-Tr)*self.BM3rdparams['dG/dT'] )            
            r = self.B_M_3rd(Pressure-Pth);V=self.BMparams['V_0(cm³/mol)']/r ;
            f = 0.5 * (pow(r, 2. / 3.) - 1.0)
            K=(1+2.*f)**(5./2.)*(K_0+(3.*K_0*self.BMparams["k'"]-5.*K_0)*f+(27./2.)*(K_0*self.BMparams["k'"]-4.*K_0)*f*f)      


            
        elif order ==4:
            K_0 = self.BMparams['K(Pa)']-((Temperature-Tr)*self.BMparams['dK/dT'] )
            G_0 = self.BMparams['G(Pa)']-((Temperature-Tr)*self.BMparams['dG/dT'] )            
            r,B1,B2,B3,B4  = self.B_M_4rd(Pressure-Pth);V=self.BMparams['V_0(cm³/mol)']/r    
            f = 0.5 * (pow(r, 2. / 3.) - 1.0)    
            r = 1./r
            K=(1+2.*f)**(5./2.)*(K_0+(3.*K_0*self.BMparams["k'"]-5.*K_0)*f+(27./2.)*(K_0*self.BMparams["k'"]-4.*K_0)*f*f)        

        

        G=(1+2.*f)**(5./2.)*(G_0+(3.*K_0*self.BMparams["G'"]-5.*G_0)*f+(6.*K_0*self.BMparams["G'"]-24.*K_0-14.*G_0+4.5*K_0*self.BMparams["k'"])*f*f)#+  (Temperature-Tr)*self.BMparams["G'"]T 
        Rho =self.molar_mass/V 
        Vp=np.sqrt((K+4.*G/3.)/Rho)/1000.
        Vs=np.sqrt(G/Rho)/1000.;#print Vs
        return Vp,Vs,Rho/1000.  ,K,G               



class Olivine(EOS):
    
    def __init__(self):
        super(EOS,self).__init__()
        self.name = 'test'
        self.formula = 'Mg2Si2O4'
        
        self.n = sum(dictionarize_formula(self.formula).values())
        self.molar_mass = formula_mass(dictionarize_formula(self.formula),atomic_masses)
        
        self.HPparams = {
            #'H_0': -2172590.0,
            'S_0(J/K)': 95.1,
            'V_0(cm³/mol)': 4.366e-05,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 2.85e-05,
            'K(Pa)': 1.279555e+11,
            'G(Pa)':81599990000.0,
            "k'": 4.21796,
            "G'":1.46257,
            "K''": -3e-11
            }
        
        self.SLBparams = {
            #'F_0': -2055403.0,
            'V_0(cm³/mol)': 4.3603e-05,
            'K(Pa)': 1.279555e+11,
            "k'": 4.21796,
            'Debye_0(K)': 809.1703,
            'grueneisen_0': 0.99282,
            'q_0': 2.10672,
            'G(Pa)': 81599990000.0,
            "G'": 1.46257,
            'eta_s_0': 2.29972,
            }        

        self.BM2rdparams = {
            #'H_0': -2172590.0,
            #'S_0(J/K)': 95.1,
            'V_0(cm³/mol)': 4.366e-05,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 2.85e-05,
            'K(Pa)': 1.279555e+11,
            'G(Pa)':81599990000.0,
            "k'": 4.0,
            "G'":1.46257,
            "K''": 0,
            'dK/dT':0,
            'dG/dT':0,
            }

        self.BM3rdparams = {
            #'H_0': -2172590.0,
            #'S_0(J/K)': 95.1,
            'V_0(cm³/mol)': 4.366e-05,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 2.85e-05,
            'K(Pa)': 1.279555e+11,
            'G(Pa)':81599990000.0,
            "k'": 4.21796,
            "G'":1.46257,
            "K''": 0,
            'dK/dT':0,
            'dG/dT':0,
            }

        self.BMparams = {
            #'H_0': -2172590.0,
            #'S_0(J/K)': 95.1,
            'V_0(cm³/mol)': 4.366e-05,
            #'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0(K-1)': 2.85e-05,
            'K(Pa)': 1.279555e+11,
            'G(Pa)':81599990000.0,
            "k'": 4.21796,
            "G'":1.46257,
            "K''": 0,
            'dK/dT':0,
            'dG/dT':0,
            } 
        
        
test=Olivine()
if __name__ == "__main__":
    test=Olivine()       

    #print ( test.B_M_4rd(pressure=1e11) )
    #print (test.Stix_Compressive(P=1e11,T=298))
    #print (test.HP_Pth(1000,Tr=300))
    #print (test.Thermal_expansion_pressure(1000,Tr=300))
    #print (test.Stix_thermal_pressure(1e5,1000))
    
    print (test.Stix_Vp_Vs(1e11,1000))
    print (test.HP_Vp_Vs(1e11,1000))
    print (test.Thermal_Vp_Vs(1e11,1000,order=2))
    print (test.Thermal_Vp_Vs(1e11,1000,order=3))
    print (test.Thermal_Vp_Vs(1e11,1000,order=4))
    
#==============================================================================
#     import matplotlib.pyplot as plt
#     #import burnman
#     P = np.linspace(0,1e11,300)
#     T = np.linspace(1300,2000,300)
#     Depth = np.linspace(0,800,300)
#     
#     v_Stix = np.array([test.Stix_Vp_Vs(P[i],T[i]) for i in range(300)])
#     v_HP = np.array([test.HP_Vp_Vs(P[i],T[i]) for i in range(300)])
#     v_thermal = np.array([test.Thermal_Vp_Vs(P[i],T[i]) for i in range(300)])
#     v_thermal2 = np.array([test.Thermal_Vp_Vs(P[i],T[i],order=2) for i in range(300)])
#    # v_thermal4 = np.array([test.Thermal_Vp_Vs(P[i],T[i],order=4) for i in range(300)])
#     
#     plt.plot(Depth,v_Stix[:,1],label='Stix')
#     plt.plot(Depth,v_HP[:,1],label = 'HP')
#     plt.plot(Depth,v_thermal[:,1],label='Thermal 3rd')
#     plt.plot(Depth,v_thermal2[:,1],label='Thermal 2rd')
#     #plt.plot(Depth,v_thermal4[:,1],label='Thermal 4rd')
#     plt.legend(loc=2)
#     plt.xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
#     plt.ylabel('$V_s$ (km/s)', fontname="Times New Roman",fontsize=10)
#     plt.title('OL vs')
#     
#     
#     plt.plot(Depth,v_Stix[:,2],label='Stix')
#     plt.plot(Depth,v_HP[:,2],label = 'HP')
#     plt.plot(Depth,v_thermal[:,2],label='Thermal 3rd')
#     plt.plot(Depth,v_thermal2[:,2],label='Thermal 2rd')
#     #plt.plot(Depth,v_thermal4[:,2],label='Thermal 4rd')
#     plt.legend(loc=2)
#     plt.xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
#     plt.ylabel('Density (kg/m3)', fontname="Times New Roman",fontsize=10)
#     plt.title('OL Densiy')
# 
#     plt.plot(Depth,v_Stix[:,0],label='Stix')
#     plt.plot(Depth,v_HP[:,0],label = 'HP')
#     plt.plot(Depth,v_thermal[:,0],label='Thermal 3rd')
#     plt.plot(Depth,v_thermal2[:,0],label='Thermal 2rd')
#     #plt.plot(Depth,v_thermal4[:,0],label='Thermal 4rd')
#     plt.legend(loc=2)
#     plt.xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
#     plt.ylabel('$V_p$ (km/s)', fontname="Times New Roman",fontsize=10)
#     plt.title('OL vp')    
#==============================================================================
