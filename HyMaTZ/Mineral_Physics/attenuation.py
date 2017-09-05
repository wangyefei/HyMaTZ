# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:53:57 2017

@author: Fei
"""
import numpy as np
from scipy.constants import gas_constant
class Anelastic():
    '''
    This class calculate anelastic based on different models 
    '''
    def __init__(self):
        pass
    
    def Melting(self,Pressure):
        # change melt point based on pressure
        return 1800 + Pressure/1e4*110
    
    def Anelastic(self,P=np.array([0.5*1e5]),T=np.array([1300]),B=0.5,w=1,a=0.2,g1=20,g2=10,Vp=np.array([9.]),Vs=np.array([3.])):
        '''
        Cammarano, et al 2003
        '''
        g = np.linspace(g1,g2,len(P))
        Tm = self.Melting(Pressure=P)
        Qs=B*w**a*np.exp(a*g*(Tm/T))*np.ones(len(Vs));b=1/Qs
        L=(4./3.)*((Vs/Vp)**2)
        bp=(L*b);Qp=1/bp
        c=np.tan(np.pi*a/2)
        return Qs,Qp,1-(b/2*c),1-(bp/2*c)
    
    def Anelastic1(self,P=np.array([0.5*1e5]),T=np.array([1300]),A=1.48*1e-1,w=0.05,a=0.15,H=500,V=20,Vp=np.array([9.]),Vs=np.array([3.])):
        '''
        Goes & Vacher  (2000)
        '''
        Qs = A*(w**a) * np.exp(a*(H*1e3+P*V*1e-6*1e5)/(gas_constant*T))
        L=(4./3.)*(Vs*1.0/Vp)**2
        bp=(L/Qs);Qp=1/bp        
        c=np.tan(np.pi*a/2)
        b=1/Qs  
        return Qs,Qp,1-(b/2*c),1-(bp/2*c)#,(H-P*V*1e5)/(gas_constant*T)
    
    def Anelastic2(self,A=2.4*1e-7,B=2.4*1e-7,Coh=6000,w=1,alpha=0.1,beta=25,Tm=1500,T=1300,Vp=np.array([9.]),Vs=np.array([3.])):
        '''
        Karato & Vacher (1998)
        '''
        b = ((A+B*Coh)/w)**alpha *np.exp(-alpha*beta*Tm/T)*np.ones(len(Vs))
        Qs=1/b
        L=(4./3.)*(Vs/Vp)**2
        bp=(L/Qs);Qp=1/bp        
        c=np.tan(np.pi*alpha/2)
        return Qs,Qp,1-(b/2*c),1-(bp/2*c)

Anelastic = Anelastic()
if __name__ == "__main__":

    
    # Goes 2000 please try this after using Main.py
    a1,b1,c1,d1 = Anelastic.Anelastic1(P=GUI.sys.Pressure,T=GUI.sys.Temperature,A=1.48*1e-1,w=0.05,a=0.15,H=500,V=20,Vp=GUI.Vp,Vs=GUI.Vs)
    a2,b2,c2,d2 = Anelastic.Anelastic1(P=GUI.sys.Pressure,T=GUI.sys.Temperature,A=2.0*1e-4,w=0.05,a=0.25,H=584,V=21,Vp=GUI.Vp,Vs=GUI.Vs)
    file = open('Goes (2000).txt','w')
    file.write('D,P,T,Qs1,Qs2')
    for i in range(len(a1)):
        file.write(str(GUI.sys.Depth[i]));file.write(' ')
        file.write(str(GUI.sys.Pressure[i]));file.write(' ')
        file.write(str(GUI.sys.Temperature[i]));file.write(' ')
        file.write(str(a1[i]));file.write(' ')
        file.write(str(a2[i]));file.write(' ')
        file.write('\n')
    file.close()
    img = plt.imread('goes(2000).png')
    fig, ax = plt.subplots()
    ax.imshow(img, extent=[0, 400, 0, 300])
    ax.set_xscale("log", nonposx='clip')
    
    # Cammarano 2003
    a1,b1,c1,d1 = Anelastic.Anelastic(P=GUI.sys.Pressure,T=GUI.sys.Temperature,B=0.5,w=1,a=0.2,g1=20,g2=10,Vp=GUI.Vp,Vs=GUI.Vs)
    a2,b2,c2,d2 = Anelastic.Anelastic(P=GUI.sys.Pressure,T=GUI.sys.Temperature,B=0.8,w=1,a=0.2,g1=20,g2=10,Vp=GUI.Vp,Vs=GUI.Vs)
    a3,b3,c3,d3 = Anelastic.Anelastic(P=GUI.sys.Pressure,T=GUI.sys.Temperature,B=1.1,w=1,a=0.2,g1=20,g2=10,Vp=GUI.Vp,Vs=GUI.Vs)
    a4,b4,c4,d4 = Anelastic.Anelastic(P=GUI.sys.Pressure,T=GUI.sys.Temperature,B=0.035,w=1,a=0.2,g1=30,g2=15,Vp=GUI.Vp,Vs=GUI.Vs)
    a5,b5,c5,d5 = Anelastic.Anelastic(P=GUI.sys.Pressure,T=GUI.sys.Temperature,B=0.056,w=1,a=0.2,g1=30,g2=15,Vp=GUI.Vp,Vs=GUI.Vs)
    file = open('Cammarano (2003).txt','w')
    file.write('D,P,T,Qs1,Qs2')
    for i in range(len(a1)):
        file.write(str(GUI.sys.Depth[i]));file.write(' ')
        file.write(str(GUI.sys.Pressure[i]));file.write(' ')
        file.write(str(GUI.sys.Temperature[i]));file.write(' ')
        file.write(str(a1[i]));file.write(' ')
        file.write(str(a2[i]));file.write(' ')
        file.write(str(a3[i]));file.write(' ')
        file.write(str(a4[i]));file.write(' ')
        file.write(str(a5[i]));file.write(' ')
        file.write('\n')
    file.close() 
    
    img = plt.imread('cammarano(2003).PNG')
    fig, ax = plt.subplots()
    ax.imshow(img, extent=[0, 600, 800, 100])
    ax.plot(a1,GUI.sys.Depth)
    ax.plot(a2,GUI.sys.Depth)
    ax.plot(a3,GUI.sys.Depth)
    ax.plot(a4,GUI.sys.Depth)
    ax.plot(a5,GUI.sys.Depth)
    ax.set_ylabel('Depth (km)', fontname="Times New Roman",fontsize=10)
    ax.set_xlabel('Qs', fontname="Times New Roman",fontsize=10)
    fig.savefig('1.png',dpi=600)
   