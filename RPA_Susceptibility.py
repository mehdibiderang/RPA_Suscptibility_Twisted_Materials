#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 14:08:41 2024

@author: mehdi
"""

#==============================================================================

import numpy as np
from numpy import linalg as LA
import Interaction_Mat as IntHam

from parameters import *

#==============================================================================

def loaddata_func(file_name):
    
    out=np.loadtxt(file_name,float)
    
    return out

def Bare_Kappa_Spin_Basis_func(file_name,num_q):
    
    Bare_Data=loaddata_func(file_name)
        
    BKappa=np.empty((2,2,2,2),float)
    
    qx=Bare_Data[num_q,0]
    qy=Bare_Data[num_q,1]
    
    number=2
    
    for s1 in range(2):
        for s2 in range(2):
            for s3 in range(2):
                for s4 in range(2):                  
                    
                    BKappa[s1,s2,s3,s4]=Bare_Data[num_q,number]
                    
                    number+=1
                    
    Kappa_Bare=[[BKappa[0,0,0,0],BKappa[0,1,0,0],BKappa[0,0,1,0],BKappa[0,1,1,0]]
               ,[BKappa[0,0,0,1],BKappa[0,1,0,1],BKappa[0,0,1,1],BKappa[0,1,1,1]]
               ,[BKappa[1,0,0,0],BKappa[1,1,0,0],BKappa[1,0,1,0],BKappa[1,1,1,0]]
               ,[BKappa[1,0,0,1],BKappa[1,1,0,1],BKappa[1,0,1,1],BKappa[1,1,1,1]]]
                    
    return qx,qy,Kappa_Bare

def RPA_calculation_func(file_name,U,num_q):
    
    qx,qy,Bare_Kappa=Bare_Kappa_Spin_Basis_func(file_name,num_q)
    
    U_Mat=IntHam.U_Mat(U)

    Id=np.identity(4,int)
        
    Coeff=LA.inv(Id-np.dot(Bare_Kappa,U_Mat))
    
    Output_Kappa_RPA=np.dot(Coeff,Bare_Kappa)
    
    return qx,qy,Output_Kappa_RPA


def RPA_spin_basis_from_Matrix_2_List_func(file_name,U,num_q):
    
    Output_Spin_RPA_List=np.empty((1,10),float)
    
    qx,qy,RPA_Matrix=RPA_calculation_func(file_name,U,num_q)
    
    Output_Spin_RPA_List[0,0]=qx
    Output_Spin_RPA_List[0,1]=qy
       
    Output_Spin_RPA_List[0,2]=RPA_Matrix[0,0]    # uuuu
    Output_Spin_RPA_List[0,3]=RPA_Matrix[0,3]    # uddu
    Output_Spin_RPA_List[0,4]=RPA_Matrix[1,1]    # udud
    Output_Spin_RPA_List[0,5]=RPA_Matrix[1,2]    # uudd
    Output_Spin_RPA_List[0,6]=RPA_Matrix[2,1]    # dduu
    Output_Spin_RPA_List[0,7]=RPA_Matrix[2,2]    # dudu
    Output_Spin_RPA_List[0,8]=RPA_Matrix[3,0]    # duud
    Output_Spin_RPA_List[0,9]=RPA_Matrix[3,3]    # dddd
    
    return Output_Spin_RPA_List

    
def RPA_Spatial_ListForm_func(file_name,U,num_q):
    
    Output_Spatial_RPA_List=np.empty((1,10),complex)
    
    qx,qy,RPA_Kappa=RPA_calculation_func(file_name,U,num_q)
    
    Output_Spatial_RPA_List[0,0]=qx
    Output_Spatial_RPA_List[0,1]=qy
    
    Output_Spatial_RPA_List[0,2]=RPA_Kappa[0,0]+RPA_Kappa[3,3]    # Kappa_00
    Output_Spatial_RPA_List[0,3]=RPA_Kappa[0,0]-RPA_Kappa[3,3]    # Kappa_0z
    Output_Spatial_RPA_List[0,4]=RPA_Kappa[1,2]+RPA_Kappa[2,1]    # Kappa_xx
    Output_Spatial_RPA_List[0,5]=-1j*(RPA_Kappa[1,2]-RPA_Kappa[2,1])    # Kappa_xy
    Output_Spatial_RPA_List[0,6]=1j*(RPA_Kappa[1,2]-RPA_Kappa[2,1])    # Kappa_yx
    Output_Spatial_RPA_List[0,7]=RPA_Kappa[1,2]+RPA_Kappa[2,1]    # Kappa_yy
    Output_Spatial_RPA_List[0,8]=RPA_Kappa[0,0]-RPA_Kappa[3,3]    # Kappa_z0
    Output_Spatial_RPA_List[0,9]=RPA_Kappa[0,0]+RPA_Kappa[3,3]    # Kappa_zz
    
    return Output_Spatial_RPA_List

def RPA_Kappa_00_func(file_name,U,num_q):
    
    Output_RPA_Kappa_00=np.empty((1,3),float)
    
    RPA_Kappa=RPA_Spatial_ListForm_func(file_name,U,num_q)
    
    Output_RPA_Kappa_00[0,0]=RPA_Kappa[0,0].real
    Output_RPA_Kappa_00[0,1]=RPA_Kappa[0,1].real
    Output_RPA_Kappa_00[0,2]=RPA_Kappa[0,2].real
    
    return Output_RPA_Kappa_00

def RPA_Kappa_zz_func(file_name,U,num_q):
    
    Output_RPA_Kappa_zz=np.empty((1,3),float)
    
    RPA_Kappa=RPA_Spatial_ListForm_func(file_name,U,num_q)
    
    Output_RPA_Kappa_zz[0,0]=RPA_Kappa[0,0].real
    Output_RPA_Kappa_zz[0,1]=RPA_Kappa[0,1].real
    Output_RPA_Kappa_zz[0,2]=RPA_Kappa[0,9].real
    
    return Output_RPA_Kappa_zz

def RPA_Kappa_pm_func(file_name,U,num_q):
    
    Output_RPA_Kappa_pm=np.empty((1,3),float)
    
    RPA_Kappa=RPA_Spatial_ListForm_func(file_name,U,num_q)
    
    Output_RPA_Kappa_pm[0,0]=RPA_Kappa[0,0].real
    Output_RPA_Kappa_pm[0,1]=RPA_Kappa[0,1].real
    Output_RPA_Kappa_pm[0,2]=0.5*(RPA_Kappa[0,4]+RPA_Kappa[0,7]-1j*(RPA_Kappa[0,5]-RPA_Kappa[0,6])).real
    
    return Output_RPA_Kappa_pm
    
    
#==============================================================================   
    
     
    