# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 16:27:40 2016
@author: Fei

Update: 10/3/2016
"""




from __future__ import absolute_import
import re
import os
from fractions import Fraction
import numpy as np
from scipy.odr import RealData,Model,ODR
from scipy.interpolate import interp1d


'''
These functions are from Burnman, but not being use in HyMaTZ anymore.
'''
def read_masses():
    """
    This function read atom's mass from a txt file and stored it into a dictionary
    datafile from Burnman
    """
    directory=os.path.dirname(__file__) 
    newfile_name = os.path.join( directory , "EXPDATA","atomic_masses.txt")
    d = {}
    with open(newfile_name) as f:
        for line in f:
           (key, val) = line.split()
           d[(key)] = float(val)
    return d
    
atomic_masses = read_masses()

def dictionarize_formula(formula):
    """
    A function to read a chemical formula string and
    convert it into a dictionary
    This function is from Burnman
    """
    f = dict()
    elements = re.findall('[A-Z][^A-Z]*', formula)
    for element in elements:
        element_name = re.split('[0-9][^A-Z]*', element)[0]
        element_atoms = re.findall('[0-9][^A-Z]*', element)
        if len(element_atoms) == 0:
            element_atoms = Fraction(1.0)
        else:
            element_atoms = Fraction(element_atoms[0])
        f[element_name] = f.get(element_name, 0.0) + element_atoms

    return f


def formula_mass(formula, atomic_masses):
    """
    This function calculate the Mineral's mass
    This function is from Burnman
    """
    mass = sum(
        formula[element] * atomic_masses[element] for element in formula)
    return mass


'''
regression plane
'''

def func(beta,data):
    """
    This function used to do the linear regression
    """
    x,y = data
    try:
        a,b,c = beta.astype(float)
    except:
        a,b,c,d = beta.astype(float)
    return a*x+b*y+c 


def Regression_Plane1(x,y,z):
    '''
    function: z=a*x+b*y+c 
    z, x, y are merusement  data
    a,b,c are parameter needed to be determine
    return (a,b,c) and their standard deviation
    '''
    data=RealData([x,y],z)
    def func(beta,data):
        x,y = data
        a,b,c = beta
        return a*x+b*y+c    
    model = Model(func)
    odr = ODR(data, model,[0.0,0.0,0.0])
    res = odr.run()
    return res.beta, res.sd_beta 

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)
    
    
def Interprate(low,high,indepedent,dependent,num=100):
    '''
    interprate two serious of data  (indepedent and dependent)
    '''
    if low>= indepedent[0] and high <= indepedent[-1]:
        x = np.linspace(low,high, num, endpoint=True)
        f = interp1d(indepedent, dependent)
        y=f(x)
    else:
        x=np.zeros(num)
        x[1:-1] = np.linspace(indepedent[0], indepedent[-1], num-2, endpoint=True)
        x[0]=low;x[-1]=high
        f=interp1d(indepedent,dependent)
        y=np.zeros(num)
        y[1:-1]=f(x[1:-1])
        y[0]=y[1]-(y[2]-y[1])/(x[2]-x[1])*(x[1]-low)
        y[-1]=y[-2]+(y[-2]-y[-3])/(x[-2]-x[-3])*(high-x[-2])
        
    return x,y
    
font = {'family' : 'serif',  
    'color'  : 'blue',  
    'weight' : 'normal',  
    'size'   : 10,  
    }            


        
def Interprate_single(indepedent,dependent,depth):
    '''
    same as interprate
    '''
    try:     
        x = depth
        f = interp1d(indepedent, dependent)
        return f(x)
    except:
        raise ValueError ( 'out of bound' )    

    
if __name__ == "__main__":
    '''
    test function Interprate
    '''
    #Pressure=[[0,1,2,3,4],[0,1,2,3,10]]
    #Pressure[1]=np.exp(Pressure[0])
    #f = interp1d(Pressure[0], Pressure[1])
    #x,y=Interprate(-1, 4,Pressure[0], Pressure[1],num=1000)
    #import matplotlib.pyplot as plt
    #plt.plot(x,y)
    import os
    print (os.path.dirname(os.path.realpath(__file__)))