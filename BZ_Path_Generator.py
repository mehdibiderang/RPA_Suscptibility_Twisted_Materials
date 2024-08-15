#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:04:08 2024

@author: mehdi
"""
import numpy as np
from math import pi,sqrt

#==============================================================================
# BZ_Path

def BZ_Path_Generator_func(dim_k_path):
    
    COEFF=2*pi/3
    
    
    # \Gamma to K
    kx_path_1=np.linspace(0,1,dim_k_path)
    ky_path_1=np.linspace(0,sqrt(3),dim_k_path)
    
    # K to M
    kx_path_2=np.linspace(1,0,dim_k_path)
    ky_path_2=np.linspace(sqrt(3),sqrt(3),dim_k_path)
    
    # M to Gamma
    kx_path_3=np.linspace(0,0,dim_k_path)
    ky_path_3=np.linspace(sqrt(3),0,dim_k_path)
    
    
    KX_path=np.empty(3*dim_k_path,float)
    KY_path=np.empty(3*dim_k_path,float)
    
    KX_path[0:dim_k_path]=COEFF*kx_path_1
    KX_path[dim_k_path:2*dim_k_path]=COEFF*kx_path_2
    KX_path[2*dim_k_path:3*dim_k_path]=COEFF*kx_path_3
    
    
    KY_path[0:dim_k_path]=COEFF*ky_path_1
    KY_path[dim_k_path:2*dim_k_path]=COEFF*ky_path_2
    KY_path[2*dim_k_path:3*dim_k_path]=COEFF*ky_path_3
    
    return KX_path,KY_path

#==============================================================================
