# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 23:41:26 2016
@author: Fei
"""

from __future__ import absolute_import

import numpy as np
import os
import scipy.constants as constants
from scipy import interpolate
G=constants.gravitational_constant



def read_table(filename):  
    '''
    read data from given file
    '''
    d=os.path.join(os.path.dirname(__file__),'DATA')

    table = []
    address=os.path.join(d, filename)
    file=open(address)
    next(file)
    for line in file:
        line=line.split();
        if line[0] != '#' and line[0] != '\xef' and line[0] != '\xef\xbb\xbf#':
            numbers = [float(i) for i in line]
            table.append(np.array((numbers)))
    file.close()
        
    
    return np.array(table)

def read_table_csv(filename):  
    d=os.path.join(os.path.dirname(__file__),'DATA')
    table = []
    address=os.path.join(d, filename)
    file=open(address)
    next(file)
    for line in file:
        line=line.split(',');
        if line[0] != '#' and line[0] != '\xef' and line[0] != '\xef\xbb\xbf#':
            numbers = [float(i) for i in line]
            table.append(np.array((numbers)))
    file.close()
   
    return np.array(table)
    





class IASP91():

    """
        Reads  REF/STW05 (STW105.txt, :cite:`kustowski2008`).
        See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        table = read_table(
            "iasp91.txt") # depth, radius, v_p, v_s
        table = np.array(table)
        self.table_depth = table[1:, 0]
        self.table_radius = table[1:, 1]
        self.table_vp = table[1:, 2]
        self.table_vs = table[1:, 3]


class AK135():

    """
        Reads  AK135 (ak135.txt, :cite:`kennett1995`).
        See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        table = read_table(
            "ak135.txt") # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_depth = table[1:, 0]
        self.table_radius = table[1:, 1]
        self.table_density = table[1:, 2]
        self.table_vp = table[1:, 3]
        self.table_vs = table[1:, 4]
        self.table_QG = table[1:, 5]
        self.table_QK = table[1:, 6]


class Suzan_extra():
    
    def __init__(self,depth=None,vs=None):
        self.table_depth = depth*1000.
        self.table_vs = vs*1000.

class Suzanref():
    '''
    Reads Suzan data, please see van der Lee, S., and A. Frederiksen (2005)
    '''
    def __init__(self):
        table=read_table(
                    "Suzan_ref.txt")
        table=np.array(table)
        self.table_depth = table[1:, 0]*1000.
        self.table_vs = table[1:, 1]*1000.

class Suzannul():
    '''
    Reads Suzan data, please see van der Lee, S., and A. Frederiksen (2005)
    '''
    def __init__(self):
        table=read_table(
                    "Suzannul.txt")
        table=np.array(table)
        self.table_depth = table[1:, 0]*1000.
        self.table_vs = table[1:, 1]*1000.

class Suzan_test():
    '''
    Reads Suzan data, please see van der Lee, S., and A. Frederiksen (2005)
    '''
    def __init__(self,name=None):
        table=read_table(
                    name)
        table=np.array(table)
        self.table_depth = table[1:, 0]*1000.
        self.table_vs = table[1:, 1]*1000.

class Suzanone():
    '''
    Reads Suzan data, please see van der Lee, S., and A. Frederiksen (2005)
    '''
    def __init__(self):
        table=read_table(
                    "Suzanone.txt")
        table=np.array(table)
        self.table_depth = table[1:, 0]*1000.
        self.table_vs = table[1:, 1]*1000.

class PREM_Suzan0():
    '''
    Reads Suzan data, please see van der Lee, S., and A. Frederiksen (2005)
    '''
    def __init__(self):
        table=read_table("PREM2.txt")
        table=np.array(table)
        self.table_depth = table[1:, 0] #Unit km
        self.table_gravity=table[1:, 1]
        self.table_pressure = table[1:, 2] #Gpa

    def gravity(self, depth):
        fgravity=interpolate.interp1d(prem_Suzan0.table_depth,prem_Suzan0.table_gravity)
        return fgravity(depth)        




class PREM_Suzan1():
    
    def __init__(self):
        table=read_table("premforFei.txt")
        table=np.array(table)
        self.table_depth = table[1:, 0] #Unit km
        self.table_density = table[1:, 1]*1000. #Unit kg/m3
        self.table_vp=table[1:, 2]   #km/s
        self.table_vs = table[1:,3]   # km/s


class PREM_C():
    
    def __init__(self):
        #SeismicTable.__init__(self)
        table = read_table_csv("PEMC.csv")
        table = np.array(table)
        self.table_depth = table[1:, 0] #Unit km
        self.table_density = table[1:, 1]*1000. #Unit kg/m3
        self.table_vp=table[1:, 2]   #km/s
        self.table_vs = table[1:,3]   # km/s
PEMC = PREM_C()        
        
class PREM_Suzan_sum():

    def __init__(self):
        prem_Suzan0=PREM_Suzan0()
        prem_Suzan1=PREM_Suzan1()     
        self.table_depth=table_depth=np.linspace(3, 2891.0, num=1000)
        fpressure=interpolate.interp1d(prem_Suzan0.table_depth,prem_Suzan0.table_pressure)
        fgravity=interpolate.interp1d(prem_Suzan0.table_depth,prem_Suzan0.table_gravity)
        self.table_gravity=fgravity(table_depth)
        self.table_pressure=fpressure(table_depth)

        
        fdensity=interpolate.interp1d(prem_Suzan1.table_depth,prem_Suzan1.table_density)   
        fvp=interpolate.interp1d(prem_Suzan1.table_depth,prem_Suzan1.table_vp)        
        fvs=interpolate.interp1d(prem_Suzan1.table_depth,prem_Suzan1.table_vs)          
        self.table_density=fdensity(table_depth)
        self.table_vp=fvp(table_depth)
        self.table_vs=fvs(table_depth)
    
class PREM():

    """
    Reads  PREM (1s) (prem.txt, :cite:`dziewonski1981`).
    See also :class:`burnman.seismic.SeismicTable`.
    """

    def __init__(self):
        table = read_table("prem.txt")
                                 # radius, pressure, density, v_p, v_s
        table = np.array(table)
        self.table_depth = table[1:, 0]
        self.table_radius = table[1:, 1]
        self.table_pressure = table[1:, 2]
        self.table_density = table[1:, 3]
        self.table_vp = table[1:, 4]
        self.table_vs = table[1:, 5]
        self.table_QK = table[1:, 6]
        self.table_QG = table[1:, 7]
        self.table_gravity=table[1:, 8]        



#PREM = PREM()
AK135=AK135()
IASP91=IASP91()
#Suzan_ref=Suzanref()
#Suzan_nul=Suzannul()
#Suzan_one=Suzanone()
prem_Suzan0=PREM_Suzan0()
prem_Suzan1=PREM_Suzan1()
PREM_Suzan=PREM_Suzan_sum()
#Suzan_test1 = Suzan_test(name="text.txt")


if __name__ == "__main__": 
    pass
    #plt.plot(Suzan_nul.table_depth,Suzan_nul.table_vs)
    #plt.plot(Suzan_one.table_depth,Suzan_one.table_vs)

