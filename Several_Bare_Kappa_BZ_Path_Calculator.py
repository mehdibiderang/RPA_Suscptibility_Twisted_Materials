#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:35:21 2024

@author: mehdi
"""

#==============================================================================

import numpy as np
import Mu_Calculator as MUcalc
import Bare_Susceptibility as KappaCalc
import BZ_Path_Generator as BZGen

from math import pi

from parameters import *



#==============================================================================

def Bare_Kappa_Calculator_func(kx_BZ, ky_BZ, tt, num_filling, num_phai, Temp, dim_k, dim_k_path, dim_mu):
    
    Filling=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8]
    Phai=[0,pi/6,pi/4,pi/3,pi/2,2*pi/3,3*pi/4,5*pi/6,pi] 
    
    filling=Filling[num_filling]
    phai=Phai[num_phai]
        
    mu=MUcalc.Mu_Calculator_func(kx_BZ, ky_BZ, tt, filling, phai, Temp, dim_k, dim_mu)
    
    print("=============================")
    print("phai = ",phai)
    print("filling = ",filling)
    print("mu = ",mu)
    print("=============================")
    
    #==============================================================================
    
    dim_q=3*dim_k_path
    
    Re_Bare_Kappa=np.empty((2,2,2,2,2,2),float)
    Im_Bare_Kappa=np.empty((2,2,2,2,2,2),float)
    
    Re_Bare_Kappa_pp=np.empty((dim_q,18),float)
    Re_Bare_Kappa_pm=np.empty((dim_q,18),float)
    Re_Bare_Kappa_mp=np.empty((dim_q,18),float)
    Re_Bare_Kappa_mm=np.empty((dim_q,18),float)
    
    Im_Bare_Kappa_pp=np.empty((dim_q,18),float)
    Im_Bare_Kappa_pm=np.empty((dim_q,18),float)
    Im_Bare_Kappa_mp=np.empty((dim_q,18),float)
    Im_Bare_Kappa_mm=np.empty((dim_q,18),float)

    #==============================================================================
    
    KX_path,KY_path=BZGen.BZ_Path_Generator_func(dim_k_path)
    
    ss=[-1,1]
                    
    for ii in range(dim_q):
        
        qx=KX_path[ii]
        qy=KY_path[ii]
        
        for band in range(2):
            for bandp in range(2):
                
                for s1 in range(2):
                    for s2 in range(2):
                        for s3 in range(2):
                            for s4 in range(2):               
                
                                Bare_Kappa=KappaCalc.Bare_Kappa_func(omega, damp, qx, qy, kx_BZ, ky_BZ, tt, mu, phai, Temp, s1 , s2, s3, s4, ss[band], ss[bandp], dim_k)
                                
                                Re_Bare_Kappa[s1,s2,s3,s4,band,bandp]=Bare_Kappa.real
                                Im_Bare_Kappa[s1,s2,s3,s4,band,bandp]=Bare_Kappa.imag
        
    
        
        
        Re_Bare_Kappa_pp[ii,0]=qx
        Re_Bare_Kappa_pp[ii,1]=qy
        
        Re_Bare_Kappa_pm[ii,0]=qx
        Re_Bare_Kappa_pm[ii,1]=qy
        
        Re_Bare_Kappa_mp[ii,0]=qx
        Re_Bare_Kappa_mp[ii,1]=qy
        
        Re_Bare_Kappa_pm[ii,0]=qx
        Re_Bare_Kappa_pm[ii,1]=qy
        
        
        
        Im_Bare_Kappa_pp[ii,0]=qx
        Im_Bare_Kappa_pp[ii,1]=qy
        
        Im_Bare_Kappa_pm[ii,0]=qx
        Im_Bare_Kappa_pm[ii,1]=qy
        
        Im_Bare_Kappa_mp[ii,0]=qx
        Im_Bare_Kappa_mp[ii,1]=qy
        
        Im_Bare_Kappa_pm[ii,0]=qx
        Im_Bare_Kappa_pm[ii,1]=qy
        
        
        number=2
        
        
        for s1 in range(2):
            for s2 in range(2):
                for s3 in range(2):
                    for s4 in range(2):               
                        
                        Re_Bare_Kappa_pp[ii,number]=Re_Bare_Kappa[s1,s2,s3,s4,1,1]
                        Re_Bare_Kappa_pm[ii,number]=Re_Bare_Kappa[s1,s2,s3,s4,1,0]
                        Re_Bare_Kappa_mp[ii,number]=Re_Bare_Kappa[s1,s2,s3,s4,0,1]
                        Re_Bare_Kappa_mm[ii,number]=Re_Bare_Kappa[s1,s2,s3,s4,0,0]
                        
                        Im_Bare_Kappa_pp[ii,number]=Im_Bare_Kappa[s1,s2,s3,s4,1,1]
                        Im_Bare_Kappa_pm[ii,number]=Im_Bare_Kappa[s1,s2,s3,s4,1,0]
                        Im_Bare_Kappa_mp[ii,number]=Im_Bare_Kappa[s1,s2,s3,s4,0,1]
                        Im_Bare_Kappa_mm[ii,number]=Im_Bare_Kappa[s1,s2,s3,s4,0,0]
                        
                        number+=1
    
    #==============================================================================
    #
    np.savetxt("Re_Bare_Kappa_pp_"+str(num_filling)+"_"+str(num_phai)+".dat",Re_Bare_Kappa_pp, fmt='%.8e', delimiter=' ', newline='\n') 
    np.savetxt("Re_Bare_Kappa_pm_"+str(num_filling)+"_"+str(num_phai)+".dat",Re_Bare_Kappa_pm, fmt='%.8e', delimiter=' ', newline='\n') 
    np.savetxt("Re_Bare_Kappa_mp_"+str(num_filling)+"_"+str(num_phai)+".dat",Re_Bare_Kappa_mp, fmt='%.8e', delimiter=' ', newline='\n') 
    np.savetxt("Re_Bare_Kappa_mm_"+str(num_filling)+"_"+str(num_phai)+".dat",Re_Bare_Kappa_mm, fmt='%.8e', delimiter=' ', newline='\n') 
    #

#==============================================================================




