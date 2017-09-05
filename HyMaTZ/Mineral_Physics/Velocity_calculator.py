# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 21:13:56 2017

@author: Fei
"""
import numpy as np



class Velocity_calculator():
    '''
    This calss contain methods to calcualte velocity 
    '''
    
    def __init__(self,K_list,G_list,Rho_list,Volume_fraction_list):
        #super(Velocity_calculator,self).__init__()
        self.K=np.array(K_list)
        self.G=np.array(G_list)
        self.Volume_fracition=np.array(Volume_fraction_list)/sum(np.array(Volume_fraction_list))
        #print (sum(self.Volume_fracition))
        self.Rho=sum(np.array(Rho_list)*self.Volume_fracition)
     
    def ReturnRho(self):
        return self.Rho

    def Voigt(self):
        K=sum(self.K*self.Volume_fracition)
        G=sum(self.G*self.Volume_fracition)
        return K,G,self.Rho
    
    def Reuss(self):
        K=1./sum(self.Volume_fracition/self.K)
        G=1./sum(self.Volume_fracition/self.G)       
        return K,G,self.Rho
    
    def Voigt_Reuss_Hill(self):
        K1,G1,r = self.Voigt()
        K2,G2,r = self.Reuss()
        K=0.5*(K1+K2)
        G=0.5*(G1+G2)
        return K,G,self.Rho

    
    '''
    This function is modified from Burnman for more detailed please visit their webpage. 
    '''
    def Hashinand_Shtrikman(self):
        K1,G1,_=self.Hashinand_Shtrikman_lower()
        K2,G2,_=self.Hashinand_Shtrikman_upper()
        K = 0.5*(K1+K2)
        G = 0.5*(G1+G2) 
        return K,G,self.Rho
    
    def Hashinand_Shtrikman_upper(self):
        if len(self.K)==0 or len(self.G)==0:
            return 0,0,self.Rho
        K_n = max(self.K)
        G_n = max(self.G)
        alpha_n = -3. / (3. * K_n + 4. * G_n)
        A_n = 0
        for i in range(len(self.Volume_fracition)):
            if self.K[i] != K_n:
                A_n += self.Volume_fracition[i]/(1./(self.K[i]-K_n)- alpha_n)

        K_upper = K_n + A_n / (1. + alpha_n * A_n)
        
        beta_n = -3. * (K_n + 2. * G_n) / (5. * G_n * (3. * K_n + 4. * G_n))
        B_n = 0
        for i in range(len(self.Volume_fracition)):
            if self.G[i] != G_n:
                B_n += self.Volume_fracition[i] / (1./(2.*(self.G[i]-G_n)) - beta_n)

        G_upper = G_n + (0.5) * B_n / (1. + beta_n * B_n)
        return K_upper,G_upper,self.Rho       
        

    def Hashinand_Shtrikman_lower(self):  
        if len(self.K)==0 or len(self.G)==0:
            return 0,0,self.Rho
        K_n = min(self.K)
        G_n = min(self.G)
        alpha_n = -3. / (3. * K_n + 4. * G_n)
        A_n = 0
        for i in range(len(self.Volume_fracition)):
            if self.K[i] != K_n:
                A_n += self.Volume_fracition[i]/(1./(self.K[i]-K_n)- alpha_n)

        K_upper = K_n + A_n / (1. + alpha_n * A_n)
        
        beta_n = -3. * (K_n + 2. * G_n) / (5. * G_n * (3. * K_n + 4. * G_n))
        B_n = 0
        for i in range(len(self.Volume_fracition)):
            if self.G[i] != G_n:
                B_n += self.Volume_fracition[i] / (1./(2.*(self.G[i]-G_n)) - beta_n)

        G_upper = G_n + (0.5) * B_n / (1. + beta_n * B_n)
        return K_upper,G_upper ,self.Rho  
    
    def Volume_average(self):
        pass


        


if __name__ == "__main__":
    V=np.array([0.1,0.2,0.3,0.4])
    K=np.array([1,2,3,4])
    G=np.array([1,2,3,4])   
    Rho=np.array([1,2,3,4])   
    KGRho = Velocity_calculator(K,G,Rho,V)
    a,b,c = KGRho.Hashinand_Shtrikman()