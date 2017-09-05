# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 20:33:49 2017
This is the old EOS　file please use EIS_modify.py
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
    
    def __init__(self):
        super(EOS,self).__init__()
    
    def Change_Name_formula(self,name,formula):
        self.params['name']=name
        self.formula = formula
    
    
    '''
    TEOS from Holand and Powell (2011)
    '''
    def HP_K_abs1(self):#K0，K1,K2 -> abc
        a=(1+self.Kprime_0)/(1+self.Kprime_0+self.K_0*self.Kdprime_0)
        b=(self.Kprime_0/self.K_0)-(self.Kdprime_0/(1+self.Kprime_0))
        c=(1+self.Kprime_0+self.K_0*self.Kdprime_0)/(self.Kprime_0**2+self.Kprime_0-self.K_0*(self.Kdprime_0**2))
        return a,b,c       

    def HP_VP0(self,Pressure):#return VP0
        a,b,c=self.HP_K_abs1()
        VP0=(1-a*(1-(1+b*Pressure)**(-c)))*self.V_0
        return VP0    

    def HP_Pth(self,Temperature,Tr=300):
        a,b,c=self.HP_K_abs1()
        TD = 10636/(self.S_0/self.n+6.44)
        u=TD/Temperature;u0=TD/Tr
        segma=(u0*u0*np.exp(u0))/(np.exp(u0)-1)**2
        Pth=(float(self.a_0)*self.K_0*TD/segma)*(1/(np.exp(u)-1)-1/(np.exp(u0)-1))
        return Pth #return Pth    
    
    def HP_VPT(self,Pressure,Temperature,Tr=298.):
        a,b,c=self.HP_K_abs1()    
        Pth=self.HP_Pth(Temperature,Tr) 
        VPT=(1-a*(1-(1+b*(Pressure-Pth))**(-c)))*self.V_0
        return VPT 

    def HP_K(self,Pressure,Temperature,Tr=298.):
        a,b,c=self.HP_K_abs1()    
        Pth=self.HP_Pth(Temperature,Tr) 
        K=self.K_0*(1+b*(Pressure-Pth))*(a+(1-a)*(1+b*(Pressure-Pth))**c)  
        return K

    def HP_G(self,Pressure,Temperature,Tr=298.,return_Rho=False):
        '''
        Shear modulus from Connolly & Kerrick (2002) and Connolly (2005)
        '''

        V=self.HP_VPT(Pressure,Temperature,Tr) 
        r=self.V_0/V               
        f = 0.5 * (pow(r, 2. / 3.) - 1.0)               
        #G=self.params['G_0'] + (Temperature-Tr)*self.params['Gprime_0T'] + Pressure*self.params['Gprime_0']
        G=(1+2.*f)**(5./2.)*(self.G_0+(3.*self.K_0*self.Gprime_0-5.*self.G_0)*f+(6.*self.K_0*self.Gprime_0-24.*self.K_0-14.*self.G_0+4.5*self.K_0*self.Kprime_0)*f*f)#-ets*(EthVT)*rho        
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
        return Vp,Vs,Rho/1000.
        
        
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
        aii=6.*self.grueneisen_0; #Eq47
        aiikk=-12.*self.grueneisen_0+36.*(self.grueneisen_0**2)-18.*self.q_0*self.grueneisen_0
        
        def fun(r): #r=rho/rho0
            #K2=(-1./self.params['K_0'])*((3-self.params['Kprime_0'])*(4-self.params['Kprime_0'])+35./9.); # Angel 2000
            f=0.5*((r)**(2./3.)-1) #eq 20 2011
            Gru=(1./6.)/(1+aii*f+0.5*aiikk*f*f)*(2*f+1)*(aii+aiikk*f)  #  eq41, 44 2005               
            Debye=abs(np.sqrt(1.+aii*f+0.5*aiikk*f*f)*self.Debye_0) #eq24 2011
                       
            #use scipy to integrae Debye function 
            EthVT =(9.*self.n*R*T)  * ((Debye/T)**(-3))*integrate.quad(Debye_int_func,0.,Debye/T)[0]
            EthVT0=(9.*self.n*R*Tr) * ((Debye/Tr)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]
                        
            Pth=(EthVT-EthVT0)*(Gru*r/self.V_0)  #
            
            if self.V_0!=0.0:
                f=0.5*(r**(2./3.)-1.)
                a=9.*self.K_0;b=27.*self.K_0*(self.Kprime_0-4.);c=0#81.0*(self.params['K_0']*K2+self.params['Kprime_0']*(self.params['Kprime_0']-7.)+143./9.); # eq 21,22 2011
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
        
        Debye=abs(np.sqrt(1+aii*f+0.5*aiikk*f*f))*self.Debye_0
        Debye_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)
        
        EthVT =(9.*self.n*gas_constant*Temperature)  * ((Debye/Temperature)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Temperature)[0]
        EthVT0=(9.*self.n*gas_constant*Tr) * ((Debye/Tr)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]
        Pth=(EthVT-EthVT0)*(Gru*r/self.V_0)
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
        R=gas_constant
        #eq 25，26 2011
        Debye_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)  # from wiki       
        r,aii,aiikk= self.Stix_Compressive(P,T,Tr)
        
        f=0.5*(r**(2./3.)-1.);
        rho=r/(self.V_0)
        V=self.V_0/r
        Rho=self.molar_mass/V
        Debye=abs(np.sqrt(1+aii*f+0.5*aiikk*f*f))*self.Debye_0
        
        EthVT=(9.*self.n*R*T)  * ((Debye/T)**(-3))*integrate.quad(Debye_int_func,0.,Debye/T)[0]
        EthVT -=(9.*self.n*R*Tr)  * ((Debye/Tr)**(-3))*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]
        
     
        #CV_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)
        #CVT=(3.*self.params['n']*R*T)*((4*(T/Debye)**(3.))*integrate.quad(CV_int_func,0.,Debye/T)[0]-3*(Debye/T)/(np.exp(Debye/T)-1))
        #CVT -=(3.*self.params['n']*R*T)*((4*(Tr/Debye)**(3.))*integrate.quad(CV_int_func,0.,Debye/Tr)[0]-3*(Debye/T)/(np.exp(Debye/T)-1))            
        CVT = 3.0 * self.n* R * (4.0 * 3*integrate.quad(Debye_int_func,0.,Debye/T)[0]/(Debye/T)**(3.) - 3.0 * (Debye/T) / (np.exp(Debye/T) - 1.0))
        CVTr = 3.0 * self.n * R * (4.0 * 3*integrate.quad(Debye_int_func,0.,Debye/Tr)[0]/(Debye/Tr)**(3.) - 3.0 * (Debye/Tr) / (np.exp(Debye/Tr) - 1.0))
        
        Gru=(1./6.)/(1+aii*f+0.5*aiikk*f*f)*(2*f+1)*(aii+aiikk*f) 

        a2s=-2.*self.grueneisen_0 - 2.*self.eta_s_0  #EQ 47
        ets=-Gru-0.5/(1.+aii*f+0.5*aiikk*f*f)*((2.*f+1.)**2.)*a2s #EQ 46
        
        q=1./9.*(18.*Gru-6.-1./2./ (1.+aii*f+0.5*aiikk*f*f) * (2.*f+1.)*(2.*f+1.)*aiikk/Gru)
        K=(1+2.*f)**(5./2.)*(self.K_0+(3.*self.K_0*self.Kprime_0-5.*self.K_0)*f+(27./2.)*(self.K_0*self.Kprime_0-4.*self.K_0)*f*f)+(Gru+1-q)*Gru*rho*(EthVT)-Gru*Gru*rho*(CVT*T-CVTr*Tr)
        a0= Gru* CVT/(K*V)  
        K*=(1.+Gru*T*a0)   #KT to KS 
        G=(1+2.*f)**(5./2.)*(self.G_0+(3.*self.K_0*self.Gprime_0-5.*self.G_0)*f+(6.*self.K_0*self.Gprime_0-24.*self.K_0-14.*self.G_0+4.5*self.K_0*self.Kprime_0)*f*f)-ets*(EthVT)*rho        
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
        return Vp,Vs,Rho/1000.        
        
    
    '''
    B_M function  Ahmad & Alkammash (2012)
    '''
    def B_M_2rd(self,pressure):
        def fun(r):    
            return pressure-1.5*self.K_0*(r**(7./3.)-r**(5./3.))
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)    
        return r
    
    def B_M_3rd(self,pressure):
        def fun(r):    
            return pressure-1.5*self.K_0*(r**(7./3.)-r**(5./3.))*(1.+0.75*(r**(2./3.)-1.)*(self.Kprime_0-4.))
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)          
        return r
    
    def B_M_4rd(self,pressure):
        B1 = (self.K_0*self.Kdprime_0) +(self.Kprime_0-4)*(self.Kprime_0-5)+(59./9.)
        B2 = (3.*self.K_0*self.Kdprime_0) + (self.Kprime_0-4)*(self.Kprime_0-13)+(129./9.)
        B3 = (3.*self.K_0*self.Kdprime_0) + (self.Kprime_0-4)*(self.Kprime_0-11)+(105./9.)
        B4 = (self.K_0*self.Kdprime_0) + (self.Kprime_0-4)*(self.Kprime_0-3) + (35./9.)
        
        def fun(r):
            return pressure- (9.*self.K_0/16.)*(-B1*r**(-5./3.)+B2*r**(-7./3.)-B3*r**(-3.)+B4*r**(-11./3.))
        try:
            r=newton(fun,1.0);
        except:
            r=bisect(fun,0.1,2.1)          
        return r     
    
   
    '''
    Thermal Pressure
    '''
    def Thermal_expansion_pressure(self,Temperature,Tr=300):
        if isinstance(self.a_0, float):
            Pth = (Temperature-Tr)*self.a_0*self.K_0
        if isinstance(self.a_0, str):
            T=Temperature
            a_T = eval(self.a_0)
            T=Tr
            a_Tr = eval(self.a_0)
            Pth=Temperature*a_T*self.K_0-Tr*a_Tr*self.K_0
        return Pth
    
   
    
    
    def Thermal_Vp_Vs(self,Pressure,Temperature,Tr=300,order=3):
        Pth = self.Thermal_expansion_pressure(Temperature,Tr)
        if order==3:
            r = self.B_M_3rd(Pressure-Pth)   
        elif order ==2:
            r = self.B_M_2rd(Pressure-Pth)
        elif order ==4:
            r = self.B_M_4rd(Pressure-Pth)            
        f = 0.5 * (pow(r, 2. / 3.) - 1.0)               
        V=self.V_0/r 
        K_0 = self.K_0-((Temperature-Tr)*self.Kdtprime )
        G_0 = self.G_0-((Temperature-Tr)*self.Gdtprime )
        K=(1+2.*f)**(5./2.)*(K_0+(3.*K_0*self.Kprime_0-5.*K_0)*f+(27./2.)*(K_0*self.Kprime_0-4.*K_0)*f*f)        
        G=(1+2.*f)**(5./2.)*(G_0+(3.*K_0*self.Gprime_0-5.*G_0)*f+(6.*K_0*self.Gprime_0-24.*K_0-14.*G_0+4.5*K_0*self.Kprime_0)*f*f)#+  (Temperature-Tr)*self.Gprime_0T 
        Rho =self.molar_mass/V 
        Vp=np.sqrt((K+4.*G/3.)/Rho)/1000.
        Vs=np.sqrt(G/Rho)/1000.;#print Vs
        return Vp,Vs,Rho/1000.        
        

class User_input(EOS):
    
    def __init__(self):
        super(EOS,self).__init__()  
        self.formula = 'Mg2SiO4'
        #formula = dictionarize_formula(self.formula)
        self.name = 'olivine'
        self.V_0 =  3.9493e-05
        self.K_0 =  1.849009e+11
        self.Kprime_0 =  4.22035
        self.Kdprime_0 =  3e-11
        self.Kdtprime = -1.86*1e-11
        self.Gdtprime = -1.36*1e-11
        self.Debye_0 =  877.7094
        self.grueneisen_0 =  1.10791
        self.q_0 =  2.3914
        self.G_0 =  1.23e+11
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  2.30461
        self.n =  7
        self.molar_mass =  formula_mass(dictionarize_formula(self.formula), atomic_masses)       
        self.H_0=-2172450.0
        self.S_0=95.1
        self.a_0='2.85e-05'
        self.Cp=[233.3, 0.001494, -603800.0, -1869.7]  
   
        
test = User_input()
            

        
        
class OL_(EOS):
    
    def __init__(self):
        super(EOS,self).__init__()  
        self.name = 'Ringwoodite'
        self.formula = '(Mg,Fe)2SiO4'
        self.V_0 =  3.9493e-05
        self.K_0 =  1.849009e+11
        self.Kprime_0 =  4.22035
        self.Kdprime_0 =  3e-11
        self.Kdtprime = -1.86*1e-11
        self.Gdtprime = -1.36*1e-11
        self.Debye_0 =  877.7094
        self.grueneisen_0 =  1.10791
        self.q_0 =  2.3914
        self.G_0 =  1.23e+11
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  2.30461
        self.n =  7
        self.molar_mass =  0.14069310000000002
        self.H_0=-2172450.0
        self.S_0=95.1
        self.a_0='2.85e-05'
        self.Cp=[233.3, 0.001494, -603800.0, -1869.7]        

class WA_(EOS):
    
    def __init__(self):
        super(EOS,self).__init__()  
        self.name = 'Ringwoodite'
        self.formula = '(Mg,Fe)2SiO4'
        self.V_0 =  3.9493e-05
        self.K_0 =  1.849009e+11
        self.Kdprime_0 =  1e-11
        self.Kprime_0 =  4.22035
        self.Kdprime_0 =  3e-11
        self.Kdtprime = -1.20*1e-11
        self.Gdtprime = -1.70*1e-11
        self.Debye_0 =  877.7094
        self.grueneisen_0 =  1.10791
        self.q_0 =  2.3914
        self.G_0 =  1.23e+11
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  2.30461
        self.n =  7
        self.molar_mass =  0.14069310000000002
        self.H_0=-2172450.0
        self.S_0=95.1
        self.a_0='2.85e-05'
        self.Cp=[233.3, 0.001494, -603800.0, -1869.7]
                

class RI_(EOS):
    
    def __init__(self):
        super(EOS,self).__init__()  
        self.name = 'Ringwoodite'
        self.formula = '(Mg,Fe)2SiO4'
        self.V_0 =  3.9493e-05
        self.K_0 =  1.849009e+11
        self.Kprime_0 =  4.22035
        self.Kdprime_0 =  3e-11
        self.Kdtprime = -2.7*1e-11
        self.Gdtprime = -1.5*1e-11
        self.Debye_0 =  877.7094
        self.grueneisen_0 =  1.10791
        self.q_0 =  2.3914
        self.G_0 =  1.23e+11
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  2.30461
        self.n =  7
        self.molar_mass =  0.14069310000000002
        self.H_0=-2172450.0
        self.S_0=95.1
        self.a_0='2.85e-05'
        self.Cp=[233.3, 0.001494, -603800.0, -1869.7]
        
 

if __name__ == "__main__":
    #example showing how to use this class
    test=OL_()       

)
    print (test.HP_Pth(1000,Tr=300))
    print (test.Thermal_expansion_pressure(1000,Tr=300))
    print (test.Stix_thermal_pressure(1e5,1000))
    
    print (test.Stix_Vp_Vs(1e11,1000))
    print (test.HP_Vp_Vs(1e11,1000))
    print (test.Thermal_Vp_Vs(1e11,1000))
    
    import matplotlib.pyplot as plt
    P = np.linspace(0,1e11,300)
    T = np.linspace(1300,2000,300)
    Depth = np.linspace(0,800,300)
    
    v_Stix = np.array([test.Stix_Vp_Vs(P[i],T[i]) for i in range(300)])
    v_HP = np.array([test.HP_Vp_Vs(P[i],T[i]) for i in range(300)])
    v_thermal = np.array([test.Thermal_Vp_Vs(P[i],T[i]) for i in range(300)])
    v_thermal2 = np.array([test.Thermal_Vp_Vs(P[i],T[i],order=2) for i in range(300)])
   # v_thermal4 = np.array([test.Thermal_Vp_Vs(P[i],T[i],order=4) for i in range(300)])
    
    plt.plot(Depth,v_Stix[:,1],label='Stix')
    plt.plot(Depth,v_HP[:,1],label = 'HP')
    plt.plot(Depth,v_thermal[:,1],label='Thermal 3rd')
    plt.plot(Depth,v_thermal2[:,1],label='Thermal 2rd')
    #plt.plot(Depth,v_thermal4[:,1],label='Thermal 4rd')
    plt.legend(loc=2)
    plt.xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
    plt.ylabel('$V_s$ (km/s)', fontname="Times New Roman",fontsize=10)
    plt.title('OL vs')
    
    
    plt.plot(Depth,v_Stix[:,2],label='Stix')
    plt.plot(Depth,v_HP[:,2],label = 'HP')
    plt.plot(Depth,v_thermal[:,2],label='Thermal 3rd')
    plt.plot(Depth,v_thermal2[:,2],label='Thermal 2rd')
    #plt.plot(Depth,v_thermal4[:,2],label='Thermal 4rd')
    plt.legend(loc=2)
    plt.xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
    plt.ylabel('Density (kg/m3)', fontname="Times New Roman",fontsize=10)
    plt.title('OL Densiy')

    plt.plot(Depth,v_Stix[:,0],label='Stix')
    plt.plot(Depth,v_HP[:,0],label = 'HP')
    plt.plot(Depth,v_thermal[:,0],label='Thermal 3rd')
    plt.plot(Depth,v_thermal2[:,0],label='Thermal 2rd')
    #plt.plot(Depth,v_thermal4[:,0],label='Thermal 4rd')
    plt.legend(loc=2)
    plt.xlabel('Depth (km)', fontname="Times New Roman",fontsize=10)
    plt.ylabel('$V_p$ (km/s)', fontname="Times New Roman",fontsize=10)
    plt.title('OL vp')    