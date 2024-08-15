#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 12:41:26 2024

@author: mehdi
"""

#==============================================================================

import numpy as np

#==============================================================================

tt=1

Temp=0.06

# BZ data

BZ_Data=np.loadtxt("BZ_Data.dat",float)
dim_k=len(BZ_Data)

kx_BZ=BZ_Data[:,0]
ky_BZ=BZ_Data[:,1]

dim_mu=101

dim_fill=31
dim_phai=29

omega=0
damp=3e-2

U=4

dim_k_path=33


#==============================================================================