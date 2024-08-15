#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 23:02:07 2024

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

def U_Mat(U):
    
    out=-U*np.kron(Pauli[2],Pauli[2])
    
    return out.real

#==============================================================================