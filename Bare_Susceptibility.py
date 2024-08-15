#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:10:50 2024

@author: mehdi
"""

#==============================================================================

import numpy as np
import Normal_Hamiltonian as NH

from parameters import *

#==============================================================================

# Pauli matrices

s0=np.identity(2,int)

sx=np.array([[0,1],[1,0]],int)
sy=np.array([[0,-1j],[1j,0]],complex)
sz=np.array([[1,0],[0,-1]],int)

Pauli=(s0,sx,sy,sz)

#==============================================================================

def Rot_U_func(kx, ky, tt, phai, band):
    
    ss=[-1,1]
    
    out=0.5*(Pauli[0]+ss[band]*NH.gk_hat_func(kx, ky, tt, phai)*Pauli[3])
    
    return out

def Coherence_Coeff_func(kx, ky, tt, phai, s1, s2, s3, s4, band, bandp):
    
    Coeff_1=Rot_U_func(kx, ky, tt, phai, band)[s1][s2]
    Coeff_2=Rot_U_func(kx, ky, tt, phai, bandp)[s3][s4]
    
    out=Coeff_1*Coeff_2
    
    return out

def Coherence_MAT_func(kx_BZ, ky_BZ, tt, phai, s1, s2, s3, s4, band, bandp,dim_k):
    
    out=np.empty(dim_k,float)
    
    for i in range(dim_k):
        
        out[i]=Coherence_Coeff_func(kx_BZ[i], ky_BZ[i], tt, phai, s1, s2, s3, s4, band, bandp)
        
    return out
                

def Kernel_func(omega, damp, qx, qy, kx, ky, tt, mu, phai, Temp, band, bandp):
    
    En=NH.energy_band_func(kx,ky,tt,mu,phai,band)
    En_q=NH.energy_band_func(kx+qx,ky+qy,tt,mu,phai,bandp)
    
    FD=NH.FD_band_func(kx,ky,tt,mu,phai,Temp,band)
    FD_q=NH.FD_band_func(kx+qx,ky+qy,tt,mu,phai,Temp,bandp)
    
    Numerator=FD_q-FD
    Denominator=omega+En-En_q+1j*damp
    
    out=Numerator/Denominator
    
    return out

def Lindhard_func(omega, damp, qx, qy, kx, ky, tt, mu, phai, Temp, s1 , s2, s3, s4, band, bandp):
    
    Coherence=Coherence_Coeff_func(kx, ky, tt, phai, s1, s2, s3, s4, band, bandp)
    Kernel=Kernel_func(omega, damp, qx, qy, kx, ky, tt, mu, phai, Temp, band, bandp)
    
    out=Coherence*Kernel
    
    return out

def Lindhard_Mat_func(omega, damp, qx, qy, kx_BZ, ky_BZ, tt, mu, phai, Temp, s1 , s2, s3, s4, band, bandp, dim_k):
    
    Coherence=Coherence_MAT_func(kx_BZ, ky_BZ, tt, phai, s1, s2, s3, s4, band, bandp, dim_k)
    Kernel=Kernel_func(omega, damp, qx, qy, kx_BZ, ky_BZ, tt, mu, phai, Temp, band, bandp)
    
    out=Coherence*Kernel
    
    return out

def Bare_Kappa_func(omega, damp, qx, qy, kx_BZ, ky_BZ, tt, mu, phai, Temp, s1 , s2, s3, s4, band, bandp, dim_k):
    
    Kernel=Lindhard_Mat_func(omega, damp, qx, qy, kx_BZ, ky_BZ, tt, mu, phai, Temp, s1 , s2, s3, s4, band, bandp, dim_k)
    
    out=sum(Kernel)/dim_k
    
    return out
    
    
#============================================================================== 

