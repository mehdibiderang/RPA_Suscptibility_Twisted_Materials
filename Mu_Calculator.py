#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:20:32 2024

@author: mehdi
"""

#==============================================================================

import numpy as np
import Normal_Hamiltonian as NH

from parameters import *

#==============================================================================

def Mu_Calculator_func(kx_BZ, ky_BZ, tt, filling, phai,Temp, dim_k, dim_mu):

    Band_Min,Band_Max=NH.BandWidth_finder_func(kx_BZ, ky_BZ, tt, phai, dim_k)
    
    mu=np.linspace(Band_Min,Band_Max,dim_mu)
    
    error_mat=np.empty(dim_mu,float)
    
    for i in range(dim_mu):
        
        fill_temp=NH.filling_func(kx_BZ, ky_BZ, tt, mu[i], phai, Temp, dim_k)
        
        error_mat[i]=abs(fill_temp-filling)
        
    error_min=min(error_mat)
    
    for i in range(dim_mu):
        
        if error_mat[i]==error_min:
            
            Mu=mu[i]
            
    return Mu
    
#==============================================================================
    
        
        
        
        

