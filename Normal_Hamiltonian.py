#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:07:37 2024

@author: mehdi
"""

#==============================================================================

import numpy as np

from parameters import *

#==============================================================================

# Pauli matrices

s0=np.identity(2,int)

sx=np.array([[0,1],[1,0]],int)
sy=np.array([[0,-1j],[1j,0]],complex)
sz=np.array([[1,0],[0,-1]],int)

Pauli=(s0,sx,sy,sz)

#==============================================================================

# Functions

def epsi_func(kx,ky,tt,mu,phai):
    
    t=abs(tt)*np.cos(phai)
    
    out=-mu-2*t*(np.cos(kx)+2*np.cos(kx/2)*np.cos(np.sqrt(3)*ky/2))
    
    return out

def gk_func(kx,ky,tt,phai):
    
    alpha=abs(tt)*np.sin(phai)
    
    out=-2*alpha*(np.sin(kx)-2*np.sin(kx/2)*np.cos(np.sqrt(3)*ky/2))
    
    return out

def gk_hat_func(kx,ky,tt,phai):
    
    gk=gk_func(kx,ky,tt,phai)
    gk_norm=abs(gk)+1e-4
    
    out=gk/gk_norm
    
    return out

def Ham_func(kx,ky,tt,mu,phai):
    
    Ham=epsi_func(kx,ky,tt,mu,phai)*Pauli[0]+gk_func(kx,ky,tt,phai)*Pauli[3]
    
    return Ham

def energy_band_func(kx,ky,tt,mu,phai,band):
    
    ss=[-1,1]
    
    epsi=epsi_func(kx,ky,tt,mu,phai)
    gk=gk_func(kx,ky,tt,phai)
    
    out=epsi+ss[band]*abs(gk)
    
    return out

def FD_band_func(kx,ky,tt,mu,phai,Temp,band):
    
    ss=[-1,1]
    
    energy=energy_band_func(kx,ky,tt,mu,phai,ss[band])
    
    length=len(energy)
    
    beta_En=energy/Temp
    
    for i in range(length):
        
        if beta_En[i]>100:
            
            beta_En[i]=100
        
    FD=1/(np.exp(beta_En)+1)
    
    return FD

def energy_mat_func(kx_BZ,ky_BZ,tt,mu,phai,dim_k):
    
    ss=[-1,1]
    
    energy=np.empty((dim_k,2),float)
    
    for band in range(2):
        
        energy[:,band]=energy_band_func(kx_BZ,ky_BZ,tt,mu,phai,ss[band])
        
    return energy

def FD_mat_func(kx_BZ,ky_BZ,tt,mu,phai,Temp,dim_k):
    
    ss=[-1,1]
    
    FD=np.empty((dim_k,2),float)
    
    for band in range(2):
        
        FD[:,band]=FD_band_func(kx_BZ,ky_BZ,tt,mu,phai,Temp,ss[band])
        
    return FD

def filling_func(kx_BZ,ky_BZ,tt,mu,phai,Temp,dim_k):
    
    FD=FD_mat_func(kx_BZ,ky_BZ,tt,mu,phai,Temp,dim_k)
    
    Filling=sum(sum(FD))/dim_k
    
    return Filling

#==============================================================================


def Extreme_Finder_func(Mat,dim_k):    
   
    Temp_Max=np.empty(dim_k,float)
    Temp_Min=np.empty(dim_k,float)
    
    for i in range(dim_k):
        
        Temp_Max[i]=max(Mat[i,:])
        Temp_Min[i]=min(Mat[i,:])
        
    Final_Max=max(Temp_Max)
    Final_Min=min(Temp_Min)
    
    return Final_Min,Final_Max 
        
      

def BandWidth_finder_func(kx_BZ,ky_BZ,tt,phai,dim_k):
    
    energy_mat=energy_mat_func(kx_BZ,ky_BZ,tt,0,phai,dim_k)
        
    En_Min,En_Max=Extreme_Finder_func(energy_mat,dim_k)
    
    return En_Min,En_Max 

#==============================================================================



    