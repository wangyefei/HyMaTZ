# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 15:59:28 2016

@author: Fei
"""


from __future__ import print_function
import numpy as np
from scipy import integrate
from scipy.optimize import newton,bisect
from scipy.constants import gas_constant


try:
    from Minerals import Mineral,NAMs
    from regression import olivine,wadsleyte,ringwoodite
    
except:
    from .Minerals import Mineral,NAMs
    from .regression import olivine,wadsleyte,ringwoodite
    
import random as random

"""
Minerals from Stixrude & Lithgow-Bertelloni 2005 and references therein
Data structure from Burnman (Cottaar etal., 2014)
"""

class Stix(Mineral):
    '''
    Base class for the finite strain-Mie-Grueneiesen-Debye equation of state detailed
    in :cite:'Stixrude2005' and 'Stixrude2011'.  The equations are
    all third order in strain.
    '''
    def __init__(self):
        Mineral.__init__(self)

        self.endmembers = [1]
        self.type='Stix'
            
    
    def Return_Rho(self,pressure=None,temperautre=None):
        self.Rho_0=self.molar_mass/self.V_0
        r,aii,aiikk = self.Compressive(pressure,temperautre)
        return self.Rho_0*r

    
    def Compressive(self,P,T,Tr=300.,term='isotropic'):
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
        
    def Gibbs_energy(self,P,T,mole_fractions=1,Tr=300):
        '''
        Given pressure (P) and temperature (T),return parameters
        Units is based on input, default is SI units
        r = V0/VPT at PT      No units
        K = bulk modulus at PT     Pa
        G = Shear modulus at PT    Pa
        V = Volume at PT
        Rho = Density at PT        
        Units is based on input, default is 
        '''
        if abs(mole_fractions) <= 1e-7:
            return 0
        r,aii,aiikk = self.Compressive(P,T)
        #Debye_int_func=lambda t:(t*t*t)/(np.exp(t)-1.)  # from wiki
        R = gas_constant
        
        
        f=0.5*(r**(2./3.)-1);#rho=r/(self.V_0)
        #V=self.V_0/r
        #Rho=self.params['molar_mass']/V
        Debye=abs(np.sqrt(1+aii*f+0.5*aiikk*f*f))*self.Debye_0
        f=0.5*(r**(2./3.)-1);rho0=1./self.V_0;#rho=(1./self.params['V_0'])*r
 
        a1=9.*self.K_0;b1=27.*self.K_0*(self.Kprime_0-4.);
        c1=81.0*(self.K_0*0.001+self.Kprime_0*(self.Kprime_0-7.)+143./9.)#biikk,biikkmm,biikkmmoo
        D_func=lambda x: (np.log(1-np.exp(-x))*x*x)
        X1=Debye/T;X2=Debye/Tr
        D11=T*(X1)**(-3)*integrate.quad(D_func,0.,X1)[0];D21=Tr*(X2)**(-3)*integrate.quad(D_func,0.,X2)[0]
        D1=(D11-D21)*9*self.n*R#/1000.
        F=self.F_0+(1./(2*rho0)*a1*f*f)+(1./(6.*rho0)*b1*f*f*f)+(1./(24.*rho0)*c1*f*f*f*f)+P*self.V_0/r+D1
        try:                
            F+=T*R*self.Fe*np.log(5)/1000.
        except:
            pass
        return F/1000 #J        

    def EOS(self,P,T,Tr=300.5,term='isotropic',return_control = 'needed'):
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
        r,aii,aiikk= self.Compressive(P,T,Tr)
        
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
        if self.type != 'NAMs':
            K*=(1.+Gru*T*a0)   #KT to KS 
        G=(1+2.*f)**(5./2.)*(self.G_0+(3.*self.K_0*self.Gprime_0-5.*self.G_0)*f+(6.*self.K_0*self.Gprime_0-24.*self.K_0-14.*self.G_0+4.5*self.K_0*self.Kprime_0)*f*f)-ets*(EthVT)*rho        
        if return_control == 'needed':
            return r,K,G,V,Rho
        else:
            return r,K,G,V,Rho,CVT,Gru,a0,Debye
    
    def Anelastic(self,w=1,a=0.2,g=20,B=0.5):
        Qs=B*w**a*np.exp(a*g)
        c=np.tan(np.pi*a/2)
        b=1/Qs
        return 1-(b/2*c)
               
        
    def Vp_Vs(self,pressure,temperature,Tr=300):
        '''
        Returns unit:[km/s]`
        Returns unit:[km/s]`
        Returns unit:[kg/m3]`
        '''
        r,K,G,V,Rho=self.EOS(pressure,temperature,Tr)
        Vp=np.sqrt((K+4.*G/3.)/Rho)/1000.
        Vs=np.sqrt(G/Rho)/1000.;#print Vs
        return Vp,Vs,Rho/1000.
    
    
    def parameters(self,pressure,temperature):
        name=self.name
        formuala = self.formula
        F_0=self.F_0
        V_0=self.V_0*1e5
        K_0=self.K_0/1e9
        Kprime_0=self.Kprime_0
        Debye_0=self.Debye_0
        grueneisen_0=self.grueneisen_0
        q_0=self.q_0
        G_0=self.G_0/1e9
        Gprime_0=self.Gprime_0
        eta_s_0=self.eta_s_0
        Vp,Vs,Rho = self.Vp_Vs(pressure,temperature)        
        return [name,formuala,Vp,Vs,Rho*1000,F_0,V_0,K_0,Kprime_0,Debye_0,grueneisen_0,q_0,G_0,Gprime_0,eta_s_0]

    def parameters_random(self,pressure,temperature,distribution = 'Gaussian'):
        '''
        Change parameters within standard deviation
        random change can base on Gaussian distribution.
        '''
        if distribution == "Gaussian":
            name=self.name
            formuala = self.formula
                
            F_0=self.random_gaussian(self.F_0, self.uncertainties['err_F_0'])           
            V_0=(self.random_gaussian(self.V_0, self.uncertainties['err_V_0']))*1e5
            K_0=(self.random_gaussian(self.K_0, self.uncertainties['err_K_0']))/1e9
            Kprime_0=self.random_gaussian(self.Kprime_0, self.uncertainties['err_K_prime_0'])
            Debye_0=self.random_gaussian(self.Debye_0, self.uncertainties['err_Debye_0'])
            grueneisen_0=self.random_gaussian(self.grueneisen_0, self.uncertainties['err_grueneisen_0'])
            q_0=self.random_gaussian(self.q_0, self.uncertainties['err_q_0'])
            G_0=self.random_gaussian(self.G_0, self.uncertainties['err_G_0'])/1e9
            Gprime_0=self.random_gaussian(self.Gprime_0, self.uncertainties['err_Gprime_0'])
            eta_s_0=self.random_gaussian(self.eta_s_0, self.uncertainties['err_eta_s_0'])
            Vp,Vs,Rho = self.Vp_Vs(pressure*1e4,temperature)  
        else:
            name=self.name
            formuala = self.formula
            F_0=self.F_0+random.uniform(-1, 1)*self.uncertainties['err_F_0']
            V_0=(self.V_0+random.random()*self.uncertainties['err_V_0'])*1e5
            K_0=(self.K_0+random.random()*self.uncertainties['err_K_0'])/1e9
            Kprime_0=self.Kprime_0+random.uniform(-1, 1)*self.uncertainties['err_K_prime_0']
            Debye_0=self.Debye_0+random.uniform(-1, 1)*self.uncertainties['err_Debye_0']
            grueneisen_0=self.grueneisen_0+random.uniform(-1, 1)*self.uncertainties['err_grueneisen_0']
            q_0=self.q_0+random.uniform(-1, 1)*self.uncertainties['err_q_0']
            G_0=(self.G_0+random.uniform(-1, 1)*self.uncertainties['err_G_0'])/1e9
            Gprime_0=self.Gprime_0+random.uniform(-1, 1)*self.uncertainties['err_Gprime_0']
            eta_s_0=self.eta_s_0+random.uniform(-1, 1)*self.uncertainties['err_eta_s_0']
            Vp,Vs,Rho = self.Vp_Vs(pressure*1e4,temperature)        
        return [name,formuala,Vp,Vs,Rho*1000,F_0,V_0,K_0,Kprime_0,Debye_0,grueneisen_0,q_0,G_0,Gprime_0,eta_s_0]
    

    def random_gaussian(self,a,b):
        try:
            c = np.random.normal(a,b)
        except:
            c = a + random.uniform(-1,1)*b
        return c
        

    def volume(self):
        '''
        return volume at 1 Atm and 298K unit: cm3 mol− 1 
        '''
        return self.V_0*1e5

        
    def change_parameters(self,parms_list):
        '''
        change minerals parameters
        '''
        self.name = parms_list[0]
        #formula = parms_list[1]
        self.V_0=float(parms_list[6])/1e5
        self.K_0=float(parms_list[7])*1e9
        self.Kprime_0=float(parms_list[8])
        self.Debye_0=float(parms_list[9])
        self.grueneisen_0=float(parms_list[10])
        self.q_0=float(parms_list[11])
        self.G_0=float(parms_list[12])*1e9
        self.Gprime_0=float(parms_list[13])
        self.eta_s_0=float(parms_list[14])  
        #self.formula=dictionarize_formula(parms_list[1])
        #formula = dictionarize_formula(formula)
        #self.n=sum(formula.values())
        #self.molar_mass=formula_mass(formula, atomic_masses)
    
    def Solve_Molar_fraction(self,weight_fraction=[0,0,0,0,0,0]):
        return np.array([1])

class anorthite (Stix):

    def __init__(self):
        self.number = 3
        self.formula = 'CaAl2Si2O8'

        self.name = 'Anorthite'     
        self.F_0 =  -4014619.0
        self.V_0 =  0.00010061
        self.K_0 =  84089150000.0
        self.Kprime_0 =  4.0
        self.Debye_0 =  752.3911
        self.grueneisen_0 =  0.39241
        self.q_0 =  1.0
        self.G_0 =  39900000000.0
        self.Gprime_0 =  1.09134
        self.eta_s_0 =  1.6254
        self.n =  13
        self.molar_mass =  0.2782072
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 4000.0,
            'err_V_0': 0.0,
            'err_K_0': 5000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 2.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 1.0,
            'err_G_0': 3000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class albite (Stix):

    def __init__(self):
        self.number = 4
        self.formula = 'NaAlSi3O8'

        self.name = 'Albite'

        
        self.F_0 =  -3718799.0
        self.V_0 =  0.000100452
        self.K_0 =  59761620000.0
        self.Kprime_0 =  4.0
        self.Debye_0 =  713.7824
        self.grueneisen_0 =  0.56704
        self.q_0 =  1.0
        self.G_0 =  36000000000.0
        self.Gprime_0 =  1.3855
        self.eta_s_0 =  1.04208
        self.n =  13
        self.molar_mass =  0.262223
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 5000.0,
            'err_V_0': 0.0,
            'err_K_0': 5000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 13.0,
            'err_grueneisen_0': 0.03,
            'err_q_0': 1.0,
            'err_G_0': 5000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class spinel (Stix):

    def __init__(self):
        self.number = 5
        self.formula = 'Mg4Al8O16'


        self.name = 'Spinel'

        
        self.F_0 =  -8667568.0
        self.V_0 =  0.000159048
        self.K_0 =  1.969428e+11
        self.Kprime_0 =  5.68282
        self.Debye_0 =  842.8104
        self.grueneisen_0 =  1.02283
        self.q_0 =  2.71208
        self.G_0 =  1.085e+11
        self.Gprime_0 =  0.37303
        self.eta_s_0 =  2.66282
        self.n =  28
        self.molar_mass =  0.5690624
        self.Fe=0

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 43.76, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 32000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 33.0,
            'err_grueneisen_0': 0.04,
            'err_q_0': 0.6,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 0.6}
        Stix.__init__(self)


class hercynite (Stix):

    def __init__(self):
        self.number = 6
        self.formula = 'Fe4Al8O16'

        self.name = 'Hercynite'

        self.F_0 =  -7324009.0
        self.V_0 =  0.000163372
        self.K_0 =  2.088965e+11
        self.Kprime_0 =  5.68282
        self.Debye_0 =  763.231
        self.grueneisen_0 =  1.21719
        self.q_0 =  2.71208
        self.G_0 =  84500000000.0
        self.Gprime_0 =  0.37303
        self.eta_s_0 =  2.768
        self.n =  28
        self.molar_mass =  0.6952224
        self.Fe=4

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 97.28, 'delta_V': 0.0}]]
        
        
        self.uncertainties = {
            'err_F_0': 35000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 32.0,
            'err_grueneisen_0': 0.07,
            'err_q_0': 1.0,
            'err_G_0': 13000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class forsterite (Stix):

    def __init__(self):
        self.number = 24
        self.formula = 'Mg2SiO4'

        self.name = 'Forsterite'

        
        self.F_0 =  -2055403.0
        self.V_0 =  4.3603e-05
        self.K_0 =  1.279555e+11
        self.Kprime_0 =  4.21796
        self.Debye_0 =  809.1703
        self.grueneisen_0 =  0.99282
        self.q_0 =  2.10672
        self.G_0 =  81599990000.0
        self.Gprime_0 =  1.46257
        self.eta_s_0 =  2.29972
        self.n =  7
        self.molar_mass =  0.14069310000000002
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 1.0,
            'err_grueneisen_0': 0.03,
            'err_q_0': 0.2,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.1}
        Stix.__init__(self)


class fayalite (Stix):

    def __init__(self):
        self.number = 25
        self.formula = 'Fe2SiO4'

        self.name = 'Fayalite'

        
        self.F_0 =  -1370519.0
        self.V_0 =  4.629e-05
        self.K_0 =  1.349622e+11
        self.Kprime_0 =  4.21796
        self.Debye_0 =  618.7007
        self.grueneisen_0 =  1.06023
        self.q_0 =  3.6466
        self.G_0 =  50899990000.0
        self.Gprime_0 =  1.46257
        self.eta_s_0 =  1.02497
        self.n =  7
        self.molar_mass =  0.20377309999999998
        self.Fe=2

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 26.76, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 2.0,
            'err_grueneisen_0': 0.07,
            'err_q_0': 1.0,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 0.6}
        Stix.__init__(self)


class mg_wadsleyite (Stix):

    def __init__(self):
        self.number = 26
        self.formula = 'Mg2SiO4'

        self.name = 'Mg_Wadsleyite'

        
        self.F_0 =  -2027837.0
        self.V_0 =  4.0515e-05
        self.K_0 =  1.686948e+11
        self.Kprime_0 =  4.3229
        self.Debye_0 =  843.4973
        self.grueneisen_0 =  1.2061
        self.q_0 =  2.0188
        self.G_0 =  1.12e+11
        self.Gprime_0 =  1.44424
        self.eta_s_0 =  2.63683
        self.n =  7
        self.molar_mass =  0.14069310000000002
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 3000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 7.0,
            'err_grueneisen_0': 0.09,
            'err_q_0': 1.0,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.2,
            'err_eta_s_0': 0.4}
        Stix.__init__(self)


class fe_wadsleyite (Stix):

    def __init__(self):
        self.number = 27
        self.formula = 'Fe2SiO4'

        self.name = 'Fe_Wadsleyite'

        
        self.F_0 =  -1364668.0
        self.V_0 =  4.28e-05
        self.K_0 =  1.68591e+11
        self.Kprime_0 =  4.3229
        self.Debye_0 =  665.4492
        self.grueneisen_0 =  1.2061
        self.q_0 =  2.0188
        self.G_0 =  72000000000.0
        self.Gprime_0 =  1.44424
        self.eta_s_0 =  1.04017
        self.n =  7
        self.molar_mass =  0.20377309999999998
        self.Fe=2

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 26.76, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 7000.0,
            'err_V_0': 0.0,
            'err_K_0': 13000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 21.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 12000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class mg_ringwoodite (Stix):

    def __init__(self):
        self.number = 28
        self.formula = 'Mg2SiO4'

        self.name = 'Mg_Ringwoodite'

        
        self.F_0 =  -2017557.0
        self.V_0 =  3.9493e-05
        self.K_0 =  1.849009e+11
        self.Kprime_0 =  4.22035
        self.Debye_0 =  877.7094
        self.grueneisen_0 =  1.10791
        self.q_0 =  2.3914
        self.G_0 =  1.23e+11
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  2.30461
        self.n =  7
        self.molar_mass =  0.14069310000000002
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 8.0,
            'err_grueneisen_0': 0.1,
            'err_q_0': 0.4,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.5}
        Stix.__init__(self)


class fe_ringwoodite (Stix):

    def __init__(self):
        self.number = 29
        self.formula = 'Fe2SiO4'

        self.name = 'Fe_Ringwoodite'

        
        self.F_0 =  -1362772.0
        self.V_0 =  4.186e-05
        self.K_0 =  2.13412e+11
        self.Kprime_0 =  4.22035
        self.Debye_0 =  677.7177
        self.grueneisen_0 =  1.27193
        self.q_0 =  2.3914
        self.G_0 =  92000000000.0
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  1.77249
        self.n =  7
        self.molar_mass =  0.20377309999999998
        self.Fe=2

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 26.76, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 7000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 8.0,
            'err_grueneisen_0': 0.23,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class enstatite (Stix):

    def __init__(self):
        self.number = 7
        self.formula = 'Mg2Si2O6'

        self.name = 'Enstatite'

        
        self.F_0 =  -2913596.0
        self.V_0 =  6.2676e-05
        self.K_0 =  1.070768e+11
        self.Kprime_0 =  7.02751
        self.Debye_0 =  812.1848
        self.grueneisen_0 =  0.78479
        self.q_0 =  3.43846
        self.G_0 =  76800000000.0
        self.Gprime_0 =  1.54596
        self.eta_s_0 =  2.50453
        self.n =  10
        self.molar_mass =  0.2007774
        self.Fe=0
        
        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.4,
            'err_Debye_0': 4.0,
            'err_grueneisen_0': 0.04,
            'err_q_0': 0.4,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.1}
        Stix.__init__(self)


class ferrosilite (Stix):

    def __init__(self):
        self.number = 8
        self.formula = 'Fe2Si2O6'

        self.name = 'Ferrosilite'

        
        self.F_0 =  -2225718.0
        self.V_0 =  6.5941e-05
        self.K_0 =  1.005386e+11
        self.Kprime_0 =  7.02751
        self.Debye_0 =  674.4769
        self.grueneisen_0 =  0.71889
        self.q_0 =  3.43846
        self.G_0 =  52000000000.0
        self.Gprime_0 =  1.54596
        self.eta_s_0 =  1.07706
        self.n =  10
        self.molar_mass =  0.2638574
        self.Fe=2

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 26.76, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 4000.0,
            'err_V_0': 0.0,
            'err_K_0': 4000000000.0,
            'err_K_prime_0': 0.5,
            'err_Debye_0': 10.0,
            'err_grueneisen_0': 0.08,
            'err_q_0': 1.0,
            'err_G_0': 5000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class mg_tschermaks (Stix):

    def __init__(self):
        self.number = 9
        self.formula = 'MgAl2SiO6'

        self.name = 'Mg_Tschermaks'

        
        self.F_0 =  -3002470.0
        self.V_0 =  5.914e-05
        self.K_0 =  1.070768e+11
        self.Kprime_0 =  7.02751
        self.Debye_0 =  783.8404
        self.grueneisen_0 =  0.78479
        self.q_0 =  3.43846
        self.G_0 =  95950860000.0
        self.Gprime_0 =  1.54596
        self.eta_s_0 =  2.49099
        self.n =  10
        self.molar_mass =  0.20234990000000003
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 9000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 24.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class ortho_diopside (Stix):

    def __init__(self):
        self.number = 10
        self.formula = 'CaMgSi2O6'

        self.name = 'Ortho_Diopside'

        
        self.F_0 =  -3015827.0
        self.V_0 =  6.8054e-05
        self.K_0 =  1.070768e+11
        self.Kprime_0 =  7.02751
        self.Debye_0 =  744.6988
        self.grueneisen_0 =  0.78479
        self.q_0 =  3.43846
        self.G_0 =  58458950000.0
        self.Gprime_0 =  1.54596
        self.eta_s_0 =  1.36161
        self.n =  10
        self.molar_mass =  0.2165504
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 3000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 9.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class diopside (Stix):

    def __init__(self):
        self.number = 13
        self.formula = 'CaMgSi2O6'


        self.name = 'Diopside'

        
        self.F_0 =  -3029531.0
        self.V_0 =  6.6039e-05
        self.K_0 =  1.122413e+11
        self.Kprime_0 =  5.23885
        self.Debye_0 =  781.6146
        self.grueneisen_0 =  0.95873
        self.q_0 =  1.52852
        self.G_0 =  67000000000.0
        self.Gprime_0 =  1.37293
        self.eta_s_0 =  1.57351
        self.n =  10
        self.molar_mass =  0.2165504
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 5000000000.0,
            'err_K_prime_0': 1.8,
            'err_Debye_0': 3.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 2.0,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class hedenbergite (Stix):

    def __init__(self):
        self.number = 14
        self.formula = 'CaFeSi2O6'


        self.name = 'Hedenbergite'

        
        self.F_0 =  -2677330.0
        self.V_0 =  6.7867e-05
        self.K_0 =  1.192555e+11
        self.Kprime_0 =  5.23885
        self.Debye_0 =  701.5851
        self.grueneisen_0 =  0.93516
        self.q_0 =  1.52852
        self.G_0 =  61000000000.0
        self.Gprime_0 =  1.17647
        self.eta_s_0 =  1.5703
        self.n =  10
        self.molar_mass =  0.24809040000000002
        self.Fe=1

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 13.38, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 45000.0,
            'err_V_0': 0.0,
            'err_K_0': 4000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 2.0,
            'err_grueneisen_0': 0.06,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class clinoenstatite (Stix):

    def __init__(self):
        self.number = 15
        self.formula = 'Mg2Si2O6'

        self.name = 'Clinoenstatite'
        
        self.F_0 =  -2905918.0
        self.V_0 =  6.25e-05
        self.K_0 =  1.122413e+11
        self.Kprime_0 =  5.23885
        self.Debye_0 =  805.0547
        self.grueneisen_0 =  0.95873
        self.q_0 =  1.52852
        self.G_0 =  79496860000.0
        self.Gprime_0 =  1.62901
        self.eta_s_0 =  1.69074
        self.n =  10.0
        self.molar_mass =  0.2007774
        self.Fe=0
        
        self.uncertainties = {
            'err_F_0': 3000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 10.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class ca_tschermaks (Stix):

    def __init__(self):
        self.number = 16
        self.formula = 'CaAl2SiO6'

        self.name = 'Ca_Tschermaks'

        
        self.F_0 =  -3120253.0
        self.V_0 =  6.3574e-05
        self.K_0 =  1.122413e+11
        self.Kprime_0 =  5.23885
        self.Debye_0 =  803.6626
        self.grueneisen_0 =  0.78126
        self.q_0 =  1.52852
        self.G_0 =  75160660000.0
        self.Gprime_0 =  1.54016
        self.eta_s_0 =  1.9672
        self.n =  10.0
        self.molar_mass =  0.2181229
        self.Fe=0

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 11.525, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 5000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 5.0,
            'err_grueneisen_0': 0.0,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class jadeite (Stix):

    def __init__(self):
        self.number = 17
        self.formula = 'NaAlSi2O6'


        self.name = 'Jadeite'

        
        self.F_0 =  -2855192.0
        self.V_0 =  6.0508e-05
        self.K_0 =  1.422873e+11
        self.Kprime_0 =  5.23885
        self.Debye_0 =  820.7623
        self.grueneisen_0 =  0.903
        self.q_0 =  0.39234
        self.G_0 =  85000000000.0
        self.Gprime_0 =  1.37398
        self.eta_s_0 =  2.18453
        self.n =  10.0
        self.molar_mass =  0.2021387
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 3000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 12.0,
            'err_grueneisen_0': 0.08,
            'err_q_0': 1.4,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class hp_clinoenstatite (Stix):

    def __init__(self):
        self.number = 11
        self.formula = 'Mg2Si2O6'


        self.name = 'HP_Clinoenstatite'

        
        self.F_0 =  -2905788.0
        self.V_0 =  6.076e-05
        self.K_0 =  1.160254e+11
        self.Kprime_0 =  6.23685
        self.Debye_0 =  824.4439
        self.grueneisen_0 =  1.12473
        self.q_0 =  0.20401
        self.G_0 =  87927170000.0
        self.Gprime_0 =  1.84119
        self.eta_s_0 =  2.14181
        self.n =  10.0
        self.molar_mass =  0.2007774
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 3000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.3,
            'err_Debye_0': 7.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 0.5,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.5}
        Stix.__init__(self)


class hp_clinoferrosilite (Stix):

    def __init__(self):
        self.number = 12
        self.formula = 'Fe2Si2O6'


        self.name = 'HP_Clinoferrosilite'

        
        self.F_0 =  -2222183.0
        self.V_0 =  6.385413e-05
        self.K_0 =  1.160254e+11
        self.Kprime_0 =  6.23685
        self.Debye_0 =  691.564
        self.grueneisen_0 =  1.12473
        self.q_0 =  0.20401
        self.G_0 =  70623090000.0
        self.Gprime_0 =  1.84119
        self.eta_s_0 =  0.79216
        self.n =  10.0
        self.molar_mass =  0.2638574
        self.Fe=2

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 26.76, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 4000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 11.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class ca_perovskite (Stix):

    def __init__(self):
        self.number = 23
        self.formula = 'CaSiO3'


        self.name = 'Ca_Perovskite'
        self.F_0 =  -1463358.0
        self.V_0 =  2.745e-05
        self.K_0 =  2.36e+11
        self.Kprime_0 =  3.9
        self.Debye_0 =  795.779
        self.grueneisen_0 =  1.88839
        self.q_0 =  0.89769
        self.G_0 =  1.568315e+11
        self.Gprime_0 =  2.22713
        self.eta_s_0 =  1.28818
        self.n =  5
        self.molar_mass =  0.1161617
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 8000.0,
            'err_V_0': 0.0,
            'err_K_0': 4000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 44.0,
            'err_grueneisen_0': 0.07,
            'err_q_0': 1.6,
            'err_G_0': 12000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class mg_akimotoite (Stix):

    def __init__(self):
        self.number = 30
        self.formula = 'MgSiO3'


        self.name = 'Mg_Akimotoite'

        self.F_0 =  -1410850.0
        self.V_0 =  2.6354e-05
        self.K_0 =  2.10706e+11
        self.Kprime_0 =  5.62088
        self.Debye_0 =  935.9778
        self.grueneisen_0 =  1.18984
        self.q_0 =  2.34514
        self.G_0 =  1.32e+11
        self.Gprime_0 =  1.57889
        self.eta_s_0 =  2.80782
        self.n = 5.0
        self.molar_mass =  0.1003887
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 4000000000.0,
            'err_K_prime_0': 0.8,
            'err_Debye_0': 12.0,
            'err_grueneisen_0': 0.13,
            'err_q_0': 0.8,
            'err_G_0': 8000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class fe_akimotoite (Stix):

    def __init__(self):
        self.number = 31
        self.formula = 'FeSiO3'


        self.name = 'Fe_Akimotoite'

        self.F_0 =  -1067598.0
        self.V_0 =  2.6854e-05
        self.K_0 =  2.10706e+11
        self.Kprime_0 =  5.62088
        self.Debye_0 =  887.8709
        self.grueneisen_0 =  1.18984
        self.q_0 =  2.34514
        self.G_0 =  1.523046e+11
        self.Gprime_0 =  1.57889
        self.eta_s_0 =  3.5716
        self.n =  5.0
        self.molar_mass =0.1319287
        self.Fe=1

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 13.38, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 21000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 120.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class corundum (Stix):

    def __init__(self):
        self.number = 32
        self.formula = 'AlAlO3'


        self.name = 'Corundum'

        
        self.F_0 =  -1582454.0
        self.V_0 =  2.5577e-05
        self.K_0 =  2.525457e+11
        self.Kprime_0 =  4.33728
        self.Debye_0 =  932.5696
        self.grueneisen_0 =  1.32442
        self.q_0 =  1.30316
        self.G_0 =  1.632e+11
        self.Gprime_0 =  1.64174
        self.eta_s_0 =  2.8316
        self.n =  5.0
        self.molar_mass = 0.1019612
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 5000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 3.0,
            'err_grueneisen_0': 0.04,
            'err_q_0': 0.2,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.2}
        Stix.__init__(self)


class pyrope (Stix):

    def __init__(self):
        self.number = 18
        self.formula = 'Mg3Al2Si3O12'


        self.name = 'Pyrope'

        
        self.F_0 =  -5936538.0
        self.V_0 =  0.00011308
        self.K_0 =  1.702396e+11
        self.Kprime_0 =  4.11067
        self.Debye_0 =  823.2102
        self.grueneisen_0 =  1.01424
        self.q_0 =  1.42169
        self.G_0 =  93699990000.0
        self.Gprime_0 =  1.35756
        self.eta_s_0 =  0.98186
        self.n =  20.0
        self.molar_mass =  0.4031273
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 10000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.3,
            'err_Debye_0': 4.0,
            'err_grueneisen_0': 0.06,
            'err_q_0': 0.5,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.2,
            'err_eta_s_0': 0.3}
        Stix.__init__(self)


class almandine (Stix):

    def __init__(self):
        self.number = 19
        self.formula = 'Fe3Al2Si3O12'


        self.name = 'Almandine'

        
        self.F_0 =  -4935516.0
        self.V_0 =  0.00011543
        self.K_0 =  1.738963e+11
        self.Kprime_0 =  4.91341
        self.Debye_0 =  741.356
        self.grueneisen_0 =  1.06495
        self.q_0 =  1.42169
        self.G_0 =  96000000000.0
        self.Gprime_0 =  1.40927
        self.eta_s_0 =  2.09292
        self.n =  20.0
        self.molar_mass =  0.4977473
        self.Fe=3

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 40.14, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 29000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 5.0,
            'err_grueneisen_0': 0.06,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class grossular (Stix):

    def __init__(self):
        self.number = 20
        self.formula = 'Ca3Al2Si3O12'


        self.name = 'Grossular'

        
        self.F_0 =  -6277935.0
        self.V_0 =  0.00012512
        self.K_0 =  1.670622e+11
        self.Kprime_0 =  3.91544
        self.Debye_0 =  822.743
        self.grueneisen_0 =  1.05404
        self.q_0 =  1.88887
        self.G_0 =  1.09e+11
        self.Gprime_0 =  1.16274
        self.eta_s_0 =  2.38418
        self.n =  20.0
        self.molar_mass =  0.4504463
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 11000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 2.0,
            'err_grueneisen_0': 0.06,
            'err_q_0': 0.2,
            'err_G_0': 4000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.1}
        Stix.__init__(self)


class mg_majorite (Stix):

    def __init__(self):
        self.number = 21
        self.formula = 'Mg4Si4O12'


        self.name = 'Mg_Majorite'

        
        self.F_0 =  -5691614.0
        self.V_0 =  0.000114324
        self.K_0 =  1.651183e+11
        self.Kprime_0 =  4.21183
        self.Debye_0 =  822.458
        self.grueneisen_0 =  0.97682
        self.q_0 =  1.53581
        self.G_0 =  84999990000.0
        self.Gprime_0 =  1.42969
        self.eta_s_0 =  1.0178
        self.n =  20.0
        self.molar_mass =  0.4015548
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 10000.0,
            'err_V_0': 0.0,
            'err_K_0': 3000000000.0,
            'err_K_prime_0': 0.3,
            'err_Debye_0': 4.0,
            'err_grueneisen_0': 0.07,
            'err_q_0': 0.5,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.2,
            'err_eta_s_0': 0.3}
        Stix.__init__(self)


class jd_majorite (Stix):

    def __init__(self):
        self.number = 22
        self.formula = 'Na2Al2Si4O12'


        self.name = 'Jd_Majorite'

        
        self.F_0 =  -5518542.0
        self.V_0 =  0.00011094
        self.K_0 =  1.770772e+11
        self.Kprime_0 =  4.11067
        self.Debye_0 =  895.914
        self.grueneisen_0 =  1.01424
        self.q_0 =  1.42169
        self.G_0 =  1.25e+11
        self.Gprime_0 =  1.35756
        self.eta_s_0 =  3.30517
        self.n =  20.0
        self.molar_mass =  0.4042774
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 14000.0,
            'err_V_0': 0.0,
            'err_K_0': 7000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 18.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 4000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class quartz (Stix):

    def __init__(self):
        self.number = 44
        self.formula = 'SiO2'

        self.name = 'Quartz'

        self.F_0 =  -858853.4
        self.V_0 =  2.367003e-05
        self.K_0 =  49547430000.0
        self.Kprime_0 =  4.33155
        self.Debye_0 =  816.3307
        self.grueneisen_0 =  0.0001
        self.q_0 =  1.0
        self.G_0 =  44856170000.0
        self.Gprime_0 =  0.95315
        self.eta_s_0 =  2.36469
        self.n =  3.0
        self.molar_mass =  0.0600843
        self.Fe=0

        self.property_modifiers = [
            ['landau', {'Tc_0': 847.0, 'S_D': 5.164, 'V_D': 1.222e-06}]]

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 31.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class coesite (Stix):

    def __init__(self):
        self.number = 45
        self.formula = 'SiO2'


        self.name = 'Coesite'

        
        self.F_0 =  -855068.5
        self.V_0 =  2.0657e-05
        self.K_0 =  1.135856e+11
        self.Kprime_0 =  4.0
        self.Debye_0 =  852.4267
        self.grueneisen_0 =  0.39157
        self.q_0 =  1.0
        self.G_0 =  61600010000.0
        self.Gprime_0 =  1.24734
        self.eta_s_0 =  2.39793
        self.n =  3.0
        self.molar_mass =  0.0600843
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 9.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class stishovite (Stix):

    def __init__(self):
        self.number = 46
        self.formula = 'SiO2'


        self.name = 'Stishovite'

        
        self.F_0 =  -818984.6
        self.V_0 =  1.4017e-05
        self.K_0 =  3.143352e+11
        self.Kprime_0 =  3.75122
        self.Debye_0 =  1107.824
        self.grueneisen_0 =  1.37466
        self.q_0 =  2.83517
        self.G_0 =  2.2e+11
        self.Gprime_0 =  1.93334
        self.eta_s_0 =  4.60904
        self.n =  3.0
        self.molar_mass =  0.0600843
        self.Fe=0

        self.property_modifiers = [
            ['landau', {'Tc_0': -4250.0, 'S_D': 0.012, 'V_D': 1e-09}]]

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 8000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 13.0,
            'err_grueneisen_0': 0.17,
            'err_q_0': 2.2,
            'err_G_0': 12000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class seifertite (Stix):

    def __init__(self):
        self.formula = 'SiO2'


        self.name = 'Seifertite'

        
        self.F_0 =  -794335.4
        self.V_0 =  1.367e-05
        self.K_0 =  3.275843e+11
        self.Kprime_0 =  4.01553
        self.Debye_0 =  1140.772
        self.grueneisen_0 =  1.37466
        self.q_0 =  2.83517
        self.G_0 =  2.274532e+11
        self.Gprime_0 =  1.76965
        self.eta_s_0 =  4.97108
        self.n =  3
        self.molar_mass =  0.0600843
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 16.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class mg_perovskite (Stix):

    def __init__(self):
        self.number = 33
        self.formula = 'MgSiO3'


        self.name = 'Mg_Perovskite'

        
        self.F_0 =  -1368283.0
        self.V_0 =  2.4445e-05
        self.K_0 =  2.505264e+11
        self.Kprime_0 =  4.14
        self.Debye_0 =  905.9412
        self.grueneisen_0 =  1.56508
        self.q_0 =  1.10945
        self.G_0 =  1.729e+11
        self.Gprime_0 =  1.69037
        self.eta_s_0 =  2.56536
        self.n =  5.0
        self.molar_mass =  0.1003887
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 3000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 5.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 0.3,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.0,
            'err_eta_s_0': 0.3}
        Stix.__init__(self)


class fe_perovskite (Stix):

    def __init__(self):
        self.number = 34
        self.formula = 'FeSiO3'


        self.name = 'Fe_Perovskite'

        
        self.F_0 =  -1040920.0
        self.V_0 =  2.5485e-05
        self.K_0 =  2.721152e+11
        self.Kprime_0 =  4.14
        self.Debye_0 =  870.8122
        self.grueneisen_0 =  1.56508
        self.q_0 =  1.10945
        self.G_0 =  1.326849e+11
        self.Gprime_0 =  1.37485
        self.eta_s_0 =  2.29211
        self.n =  5.0
        self.molar_mass =  0.1319287
        self.Fe=1

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 13.38, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 6000.0,
            'err_V_0': 0.0,
            'err_K_0': 40000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 26.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 40000000000.0,
            'err_Gprime_0': 0.0,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class al_perovskite (Stix):

    def __init__(self):
        self.number = 35
        self.formula = 'AlAlO3'


        self.name = 'Al_perovskite'

        
        self.F_0 =  -1533878.0
        self.V_0 =  2.4944e-05
        self.K_0 =  2.582e+11
        self.Kprime_0 =  4.14
        self.Debye_0 =  886.4601
        self.grueneisen_0 =  1.56508
        self.q_0 =  1.10945
        self.G_0 =  1.713116e+11
        self.Gprime_0 =  1.49706
        self.eta_s_0 =  2.47126
        self.n =  5.0
        self.molar_mass =  0.1019612
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 0.5,
            'err_Debye_0': 7.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.5}
        Stix.__init__(self)


class mg_post_perovskite (Stix):

    def __init__(self):
        self.number = 36
        self.formula = 'MgSiO3'


        self.name = 'Mg_Post_Perovskite'

        
        self.F_0 =  -1348641.0
        self.V_0 =  2.4419e-05
        self.K_0 =  2.312e+11
        self.Kprime_0 =  4.0
        self.Debye_0 =  855.8173
        self.grueneisen_0 =  1.89155
        self.q_0 =  1.09081
        self.G_0 =  1.50167e+11
        self.Gprime_0 =  1.97874
        self.eta_s_0 =  1.16704
        self.n =  5.0
        self.molar_mass =  0.1003887
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 3000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 7.0,
            'err_grueneisen_0': 0.03,
            'err_q_0': 0.1,
            'err_G_0': 4000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.2}
        Stix.__init__(self)


class fe_post_perovskite (Stix):

    def __init__(self):
        self.number = 37
        self.formula = 'FeSiO3'


        self.name = 'Fe_Post_Perovskite'

        
        self.F_0 =  -981806.9
        self.V_0 =  2.5459e-05
        self.K_0 =  2.312e+11
        self.Kprime_0 =  4.0
        self.Debye_0 =  781.3465
        self.grueneisen_0 =  1.89155
        self.q_0 =  1.09081
        self.G_0 =  1.295e+11
        self.Gprime_0 =  1.44675
        self.eta_s_0 =  1.36382
        self.n =  5.0
        self.molar_mass =  0.1319287
        self.Fe=1

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 13.38, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 21000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 52.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 5000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class al_post_perovskite (Stix):

    def __init__(self):
        self.number = 38
        self.formula = 'AlAlO3'


        self.name = 'Al_Post_Perovskite'

        
        self.F_0 =  -1377582.0
        self.V_0 =  2.3847e-05
        self.K_0 =  2.49e+11
        self.Kprime_0 =  4.0
        self.Debye_0 =  762.1951
        self.grueneisen_0 =  1.64573
        self.q_0 =  1.09081
        self.G_0 =  91965310000.0
        self.Gprime_0 =  1.81603
        self.eta_s_0 =  2.83762
        self.n = 5.0
        self.molar_mass =  0.1019612
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 4000.0,
            'err_V_0': 0.0,
            'err_K_0': 20000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 9.0,
            'err_grueneisen_0': 0.02,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.2}
        Stix.__init__(self)


class periclase (Stix):

    def __init__(self):
        self.number = 42
        self.formula = 'MgO'


        self.name = 'Periclase'

        
        self.F_0 =  -569444.6
        self.V_0 =  1.1244e-05
        self.K_0 =  1.613836e+11
        self.Kprime_0 =  3.84045
        self.Debye_0 =  767.0977
        self.grueneisen_0 =  1.36127
        self.q_0 =  1.7217
        self.G_0 =  1.309e+11
        self.Gprime_0 =  2.1438
        self.eta_s_0 =  2.81765
        self.n =  2.0
        self.molar_mass =  0.040304400000000004
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 0.0,
            'err_V_0': 0.0,
            'err_K_0': 3000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 9.0,
            'err_grueneisen_0': 0.05,
            'err_q_0': 0.2,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.2}
        Stix.__init__(self)


class wuestite (Stix):

    def __init__(self):
        self.number = 43
        self.formula = 'FeO'


        self.name = 'Wuestite'

        
        self.F_0 =  -242146.0
        self.V_0 =  1.2264e-05
        self.K_0 =  1.794442e+11
        self.Kprime_0 =  4.9376
        self.Debye_0 =  454.1592
        self.grueneisen_0 =  1.53047
        self.q_0 =  1.7217
        self.G_0 =  59000000000.0
        self.Gprime_0 =  1.44673
        self.eta_s_0 =  -0.05731
        self.n =  2.0
        self.molar_mass =  0.0718444
        self.Fe=1

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 13.38, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 1000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 21.0,
            'err_grueneisen_0': 0.13,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class mg_ca_ferrite (Stix):

    def __init__(self):
        self.number = 39
        self.formula = 'MgAl2O4'

        self.name = 'Mg_Ca_Ferrite'

        
        self.F_0 =  -2122169.0
        self.V_0 =  3.6177e-05
        self.K_0 =  2.106663e+11
        self.Kprime_0 =  4.0528
        self.Debye_0 =  838.6291
        self.grueneisen_0 =  1.31156
        self.q_0 =  1.0
        self.G_0 =  1.29826e+11
        self.Gprime_0 =  1.75878
        self.eta_s_0 =  2.1073
        self.n = 7.0
        self.molar_mass =  0.1422656
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 4000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 16.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class fe_ca_ferrite (Stix):

    def __init__(self):
        self.number = 40
        self.formula = 'FeAl2O4'


        self.name = 'Fe_Ca_Ferrite'

    
        self.F_0 =  -1790284.0
        self.V_0 =  3.7258e-05
        self.K_0 =  2.106663e+11
        self.Kprime_0 =  4.0528
        self.Debye_0 =  804.1986
        self.grueneisen_0 =  1.31156
        self.q_0 =  1.0
        self.G_0 =  1.535236e+11
        self.Gprime_0 =  1.75878
        self.eta_s_0 =  3.0268
        self.n =  7.0
        self.molar_mass =  0.1738056
        self.Fe=1

        self.property_modifiers = [
            ['linear', {'delta_E': 0.0, 'delta_S': 13.38, 'delta_V': 0.0}]]

        self.uncertainties = {
            'err_F_0': 25000.0,
            'err_V_0': 0.0,
            'err_K_0': 10000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 69.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class na_ca_ferrite (Stix):

    def __init__(self):
        self.number = 41
        self.formula = 'NaAlSiO4'


        self.name = 'Na_Ca_Ferrite'

        
        self.F_0 =  -1844129.0
        self.V_0 =  3.627e-05
        self.K_0 =  1.613385e+11
        self.Kprime_0 =  4.32479
        self.Debye_0 =  812.4769
        self.grueneisen_0 =  0.69428
        self.q_0 =  1.0
        self.G_0 =  1.220049e+11
        self.Gprime_0 =  2.07687
        self.eta_s_0 =  2.79016
        self.n = 7.0
        self.molar_mass =  0.1420544
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 11000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.1,
            'err_Debye_0': 51.0,
            'err_grueneisen_0': 0.3,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class kyanite (Stix):

    def __init__(self):
        self.number = 48
        self.formula = 'Al2SiO5'


        self.name = 'Kyanite'

        
        self.F_0 =  -2446058.0
        self.V_0 =  4.4227e-05
        self.K_0 =  1.6e+11
        self.Kprime_0 =  4.0
        self.Debye_0 =  943.1665
        self.grueneisen_0 =  0.9255
        self.q_0 =  1.0
        self.G_0 =  1.204033e+11
        self.Gprime_0 =  1.7308
        self.eta_s_0 =  2.96665
        self.n =  8.0
        self.molar_mass =  0.1620455
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 4000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 0.0,
            'err_Debye_0': 8.0,
            'err_grueneisen_0': 0.07,
            'err_q_0': 1.0,
            'err_G_0': 10000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)


class nepheline (Stix):

    def __init__(self):
        self.number = 49
        self.formula = 'NaAlSiO4'


        self.name = 'Nepheline'

        
        self.F_0 =  -1992104.0
        self.V_0 =  5.46684e-05
        self.K_0 =  53077990000.0
        self.Kprime_0 =  4.0
        self.Debye_0 =  700.9422
        self.grueneisen_0 =  0.69428
        self.q_0 =  1.0
        self.G_0 =  30700000000.0
        self.Gprime_0 =  1.33031
        self.eta_s_0 =  0.6291
        self.n =  7.0
        self.molar_mass =  0.1420544
        self.Fe=0

        self.uncertainties = {
            'err_F_0': 3000.0,
            'err_V_0': 0.0,
            'err_K_0': 1000000000.0,
            'err_K_prime_0': 1.0,
            'err_Debye_0': 13.0,
            'err_grueneisen_0': 0.03,
            'err_q_0': 1.0,
            'err_G_0': 1000000000.0,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0}
        Stix.__init__(self)
        

class input_mineral (Stix):

    def __init__(self):
        self.number = 49
        self.formula = 'NaAlSiO4'


        self.name = 'Nepheline'

        
        self.F_0 =  -1992104.0
        self.V_0 =  5.46684e-05
        self.K_0 =  53077990000.0
        self.Kprime_0 =  4.0
        self.Debye_0 =  700.9422
        self.grueneisen_0 =  0.69428
        self.q_0 =  1.0
        self.G_0 =  30700000000.0
        self.Gprime_0 =  1.33031
        self.eta_s_0 =  0.6291
        self.n =  7.0
        self.molar_mass =  0.1420544
        self.Fe=0


        Stix.__init__(self)



class OL_(Stix,NAMs):

    def __init__(self,regression=None):
        self.name = "olivine"
        self.type='NAMs'
        self.formula = 'Mg2SiO4'

        self.a=regression.Return()  

        self.name = 'Olivine'
        self.formula = '(Mg,Fe)2SiO4'
        self.F_0 =  -2055403.0
        self.V_0 =  4.3603e-05
        self.K_0 =  1.279555e+11
        self.Kprime_0 =  4.21796
        self.Debye_0 =  809.1703
        self.grueneisen_0 =  0.99282
        self.q_0 =  2.10672
        self.G_0 =  81599990000.0
        self.Gprime_0 =  1.46257
        self.eta_s_0 =  2.29972
        self.n =  7
        self.molar_mass =  0.14069310000000002
        
        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 1.0,
            'err_grueneisen_0': 0.03,
            'err_q_0': 0.2,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.1}
        if olivine.pressure_deri is not None:
            self.Kprime_0=olivine.pressure_deri[0]
            self.Gprime_0=olivine.pressure_deri[1]            
        Mineral.__init__(self)
        NAMs.__init__(self)

    def change_prime(self,kprime=1.21):
        self.Kprime_0 = kprime

        

class WA_(Stix,NAMs):

    def __init__(self,regression=None):
        self.name="wadsleyte"
        self.type='NAMs'
        self.formula = 'Mg2SiO4'
        self.a=regression.Return()     

        self.name = 'Wadsleyite'
        self.formula = '(Mg,Fe)2SiO4'
        self.F_0 =  -2027837.0
        self.V_0 =  4.0515e-05
        self.K_0 =  1.686948e+11
        self.Kprime_0 =  4.3229
        self.Debye_0 =  843.4973
        self.grueneisen_0 =  1.2061
        self.q_0 =  2.0188
        self.G_0 =  1.12e+11
        self.Gprime_0 =  1.44424
        self.eta_s_0 =  2.63683
        self.n =  7
        self.molar_mass =  0.14069310000000002

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 3000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 7.0,
            'err_grueneisen_0': 0.09,
            'err_q_0': 1.0,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.2,
            'err_eta_s_0': 0.4}
        if wadsleyte.pressure_deri is not None:
            self.Kprime_0=wadsleyte.pressure_deri[0]
            self.Gprime_0=wadsleyte.pressure_deri[1]            
        Mineral.__init__(self)
        NAMs.__init__(self)

    def change_prime(self,kprime=1.21):
        self.Kprime_0 = kprime

class RI_(Stix,NAMs):

    def __init__(self,regression=None):
        self.name="ringwoodite"
        self.type='NAMs'
        self.formula = 'Mg2SiO4'
        self.a=regression.Return()       

        self.name = 'Ringwoodite'
        self.formula = '(Mg,Fe)2SiO4'
        self.F_0 =  -2017557.0
        self.V_0 =  3.9493e-05
        self.K_0 =  1.849009e+11
        self.Kprime_0 =  4.22035
        self.Debye_0 =  877.7094
        self.grueneisen_0 =  1.10791
        self.q_0 =  2.3914
        self.G_0 =  1.23e+11
        self.Gprime_0 =  1.35412
        self.eta_s_0 =  2.30461
        self.n =  7
        self.molar_mass =  0.14069310000000002

        self.uncertainties = {
            'err_F_0': 2000.0,
            'err_V_0': 0.0,
            'err_K_0': 2000000000.0,
            'err_K_prime_0': 0.2,
            'err_Debye_0': 8.0,
            'err_grueneisen_0': 0.1,
            'err_q_0': 0.4,
            'err_G_0': 2000000000.0,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.5}
        if ringwoodite.pressure_deri is not None:
            self.Kprime_0=ringwoodite.pressure_deri[0]
            self.Gprime_0=ringwoodite.pressure_deri[1]
        Mineral.__init__(self)
        NAMs.__init__(self)
        
    def change_prime(self,kprime=1.21):
        self.Kprime_0 = kprime        

  
#==============================================================================
# class test():
#     
#     def __init__(self):
#         pass
#     
#     def EOS(self,P=None,T=None):
#         return 0,0,0,0,0
#     
# 
# 
# 
#    
# test=test()
#==============================================================================
'''
NAMs class
'''
OL=OL_(olivine)
WA=WA_(wadsleyte)
RI=RI_(ringwoodite) 
        



'''
Stix aliases
'''

ab = albite()
an = anorthite()

# LP Spinels
sp = spinel()
hc = hercynite()

# Olivine polymorphs
fo = forsterite()
fa = fayalite()
mgwa = mg_wadsleyite()
fewa = fe_wadsleyite()
mgri = mg_ringwoodite()
feri = fe_ringwoodite()

# Orthopyroxenes
en = enstatite()
fs = ferrosilite()
mgts = mg_tschermaks()
odi = ortho_diopside()

# Clinopyroxenes
di = diopside()
he = hedenbergite()
cen = clinoenstatite()
cats = ca_tschermaks()
jd = jadeite()
mgc2 = hp_clinoenstatite()
fec2 = hp_clinoferrosilite()
hpcen = hp_clinoenstatite()
hpcfs = hp_clinoferrosilite()

# Perovskites
mgpv = mg_perovskite()
mg_bridgmanite = mg_perovskite()
fepv = fe_perovskite()
fe_bridgmanite = fe_perovskite()
alpv = al_perovskite()
capv = ca_perovskite()

# Ilmenite group
mgil = mg_akimotoite()
feil = fe_akimotoite()
co = corundum()

# Garnet group
py = pyrope()
al = almandine()
gr = grossular()
mgmj = mg_majorite()
jdmj = jd_majorite()

# Quartz polymorphs
qtz = quartz()
q = quartz()
coes = coesite()
coe = coesite()
st = stishovite()
seif = seifertite()

# Post perovskites
mppv = mg_post_perovskite()
fppv = fe_post_perovskite()
appv = al_post_perovskite()

# Magnesiowuestite
pe = periclase()
wu = wuestite()

# Calcium ferrite structured phases
mgcf = mg_ca_ferrite()
fecf = fe_ca_ferrite()
nacf = na_ca_ferrite()

# Al2SiO5 polymorphs
ky = kyanite()

# Nepheline group
neph = nepheline()


if __name__ == "__main__":
    P=5.9e9;T=1673
    print (feri.Gibbs_energy(P,T))
    print (fewa.Gibbs_energy(P,T))    
    print (fa.Gibbs_energy(P,T))
        
