# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 15:35:25 2017
This file use to transform  Perple_X output data to the format HyMaTZ used 
@author: Fei
"""
import os
import numpy as np
try:
    from Stix2011data import (OL,WA,RI,ab,an,sp,hc,fo,fa,mgwa,fewa,mgri,feri,en,fs,mgts,odi,di\
                              ,he,cen,cats,jd,hpcen,hpcfs,mgpv,fepv\
                              ,alpv,capv,mgil,feil,co,py ,al,gr,mgmj ,jdmj,qtz,coes,st,seif,mppv,fppv,appv,pe,wu,mgcf \
                              ,fecf ,nacf,ky,neph,q )
    
    
    from Solidsolution import (c2c,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,ppv,Pv,Ring,Sp,Wad,Gt_maj)
except:
    from .Stix2011data import (OL,WA,RI,ab,an,sp,hc,fo,fa,mgwa,fewa,mgri,feri,en,fs,mgts,odi,di\
                          ,he,cen,cats,jd,hpcen,hpcfs,mgpv,fepv\
                          ,alpv,capv,mgil,feil,co,py ,al,gr,mgmj ,jdmj,qtz,coes,st,seif,mppv,fppv,appv,pe,wu,mgcf \
                          ,fecf ,nacf,ky,neph,q )


    from .Solidsolution import (c2c,CF,Cpx,Gt,Aki,Wus,O,Opx,Pl,ppv,Pv,Ring,Sp,Wad,Gt_maj)

'''
if anything goes wrong, check if in Stix2011data.py, every function solve correctly, and place in wright place.
'''


filename = 'E:\\perplex\\Harzburgite82'
address_input=filename+'_1.phm'
address_output2 = filename+'_2.txt'
address_output3 = filename+'_3.txt'

def Extrac_data(address_input,address_output2,address_output3):
    
#==============================================================================
#     aa=os.path.dirname(os.path.abspath('__file__'))
#     f21 = 'Input';f22='Input_2.txt'
#     f31 = 'Input';f32='Input_3.txt' 
#     f2=os.path.join( aa , 'Models',f21,f22)
#     f3=os.path.join( aa , 'Models',f31,f32)
#     print (aa)
#     print (f2)    
#==============================================================================
    
    c=np.zeros((60000,55))
    header =['Name','Counter','P(bar)','T(K)','V,J/bar/mol','H,J/mol','Gruneisen_T','Ks,bar','Gs,bar','v0,km/s','vp,km/s','vs,km/s','vp/vs','rho,kg/m3','G,J/mol','cp,J/K/mol','alpha,1/K',\
          'beta,1/bar','S,J/K/mol','n,mol','N,g','Ks_T,bar/K','Gs_T,bar/K','Ks_P','Gs_P','v0_T','vp_T','vs_T','v0_P','vp_P','vs_P','cp/cv','vol,%','wt,%','mol,%','CAO,wt%','AL2O3,wt%','NA2O,wt%',\
           'MGO,wt%','FEO,wt%','SIO2,wt%','mu_CAO,J/mol','mu_AL2O3,J/mol', 'mu_NA2O,J/mol',  'mu_MGO,J/mol','mu_FEO,J/mol','mu_SIO2,J/mol']
    
    '''
            Pi_num=0;Ti_num=2;an_num=3;ab_num=4;sp_num=5;hc_num=6;
            en_num=7;fs_num=8;mgts_num=9;odi_num=10;
            hpcen_num=11;hpcfs_num=12;
            di_num=13;he_num=14;cen_num=15;cats_num=16;
            jd_num=17;py_num=18;al_num=19;gr_num=20; mgmj_num=21;jdmj_num=22;capv_num=23;
            fo_num=24;fa_num=25;mgwa_num=26;fewa_num=27;mgri_num=28;feri_num=29;#OL=50;WA=51;RI=52
            mgil_num=30;feil_num=31;co_num=32;mgpv_num=33;fepv_num=34;alpv_num=35;mppv_num=36;fppv_num=37;appv_num=38;mgcf_num=39;fecf_num=40;
            nacf_num=41;pe_num=42;wu_num=43;qtz_num=44;coes_num=45;st_num=46;apbo_num=47;ky_num=48;neph_num=49; #apbo neph did not found
        
    '''
    
    '''
    
    'P(bar)2','T(K)3','V,J/bar/mol4','H,J/mol5',''Ks,bar7','Gs,bar8',
    'vp,km/s10','vs,km/s11','rho,kg/m312','cp,J/K/mol15','alpha,1/K16',\
    'beta,1/bar','S,J/K/mol17','n,mol18','N,g','Ks_T,bar/K','Gs_T,bar/K','Ks_P','Gs_P','v0_T','vp_T','vs_T','v0_P','vp_P',
    'vs_P','cp/cv','vol,%'32,'wt,%','mol,%',
    'CAO,wt%35','AL2O3,wt%36','NA2O,wt%37','MGO,wt%38','FEO,wt%39','SIO2,wt%40'
    '''
#==============================================================================
#     class ststem():
#         def __init__(self):
#             pass
#         def Solve_Molar_fraction(self,weight_fraction=[0,0,0,0,0,0]):
#             pass
#     system=ststem()
#==============================================================================
    
    print_lines=[]
    print_lines3=[]
    file=open(address_input)
    for i in range(13):
        next(file)
    for line in file:
        line = line.split()
        #break
        try:
            aa = np.array([float(line[35]),float(line[36]),float(line[37]),float(line[38]),float(line[39]),float(line[40])])
        except:
            print (line,'aa is wrong')
        if line[0] == 'system':
            print_line = np.zeros(61)
            print_line[0] = float(line[2]) #P
            print_line[2] = float(line[3]) #T
            print_line[53] = float(line[10]) #Vp
            print_line[54] = float(line[11]) #Vs
            print_line[55] = float(line[13]) #Density
            print_line[56] = float(line[18]) #S
            print_line[57] = float(line[5]) #H
            print_line[58] = float(line[16]) #alpha
            print_line[59] = float(line[15]) #Cp
            print_line[60] = float(line[4]) #V
            print_lines3.append([line[2],line[3],line[18],line[5],line[16],line[15],line[4],line[13]])
        else:    
            try:
                if line[0] == 'ca-pv':
                    b=capv.Solve_Molar_fraction(aa)*float(line[32])/100
                    print_line[capv.number] = b
                elif line[0] == 'C2/c':
                        a = c2c    
                        b = a.Solve_Molar_fraction(aa)*float(line[32])/100
                        c = a.endmembers[0][0].number;d=len(a.endmembers)+c
                        print_line[c:d]+=b              
                else:
                    try:
                        a = eval(line[0])
                        b = a.Solve_Molar_fraction(aa)*float(line[32])/100
                        c = a.endmembers[0][0].number;d=len(a.endmembers)+c
                        print_line[c:d]+=b
                    except:
                        a = eval(line[0])
                        b = a.Solve_Molar_fraction(aa)*float(line[32])/100
                        print_line[a.number] = b  
            except:
                print (line[0])
                print_line = np.zeros(61)
                print_lines.append(print_line)
                print_lines3.append(['0','0','0','0','0','0','0','0'])
        if line[0] == 'system':
            print_lines.append(print_line)
        #break
    file.close()
    file1=open(address_output3,'w')
    file1.write('P(bar) Depth(km) Temperautre(K) an ab sp hc en fs mgts odi hpcen hpcfs di he cen cats jd py al gr  mgmj jdmj capv fo fa mgwa fewa mgri feri mgil feil co mgpv fepv alpv mppv fppv appv mgcf fecf ')
    file1.write('nacf pe wu qtz coes st apbo ky neph  OL WA RI  Vp Vs')
    file1.write('\n')
    for i in print_lines:
        for num,j in enumerate(i):
            if num <=54:# or num>=55:
                file1.write(str(j))
                file1.write(' ')
        file1.write('\n')
    file1.close()
    
    file1=open(address_output2,'w')
    file1.write('P  T   S H alpha Cp V  Density')
    file1.write('\n')
    for i in print_lines3:
        for j in i:
            file1.write(j)
            file1.write(' ')
        file1.write('\n')
    file1.close()
    
    aa=os.path.dirname(os.path.abspath('__file__'))
    f21 = 'Input';f22='Input_2.txt'
    f31 = 'Input';f32='Input_3.txt' 
    f2=os.path.join( aa , 'Models',f21,f22)
    f3=os.path.join( aa , 'Models',f31,f32)
    print (aa)
    print (f2)
    try:
        os.remove(f2)
        os.remove(f3)
    except:
        pass
    os.rename(address_output2, f2 )
    os.rename(address_output3, f3 )
    #print (aa)
    
    
if __name__ == "__main__":
    
    filename = 'C:\\Users\\Fei\\Downloads\\Harzburgite82_1.phm'
    address_input=filename
    address_output2 = filename[:-5]+'2.txt'
    address_output3 = filename[:-5]+'3.txt'


    f2 = 'E:\\Data\\Input\\Input_2.txt'
    f3 = 'E:\\Data\\Input\\Input_3.txt'
    Extrac_data(address_input,address_output2,address_output3)  
    #import os
    #os.rename(address_output2, f2 )
    