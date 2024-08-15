#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:29:48 2024

@author: mehdi
"""

#==============================================================================

import time

start_time=time.time()

#==============================================================================

import numpy as np
import Several_Bare_Kappa_BZ_Path_Calculator as SevBareKapapa
from parameters import *


#==============================================================================

num_filling=9
num_phai=9

for ii in range(num_filling):
    for jj in range(num_phai):
        
        SevBareKapapa.Bare_Kappa_Calculator_func(kx_BZ, ky_BZ, tt, ii, jj, Temp, dim_k, dim_k_path, dim_mu)



#==============================================================================

end_time=time.time()

Time=end_time-start_time

print("Time = ",Time)

#==============================================================================        
