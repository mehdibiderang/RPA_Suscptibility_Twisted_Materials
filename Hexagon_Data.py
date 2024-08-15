#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:45:33 2020

@author: delsa
"""

import numpy as np
from matplotlib import path

# from shapely.geometry import Point, Polygon


from math import pi,sqrt#,cos,sin,exp

import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
from matplotlib import rc


#====================================================================
import time
start_time = time. time()

#====================================================================

dim_k=101

Coeff=2

kk=np.linspace(Coeff*-pi,Coeff*pi,dim_k)

dk=kk[1]-kk[0]

kx,ky=np.meshgrid(kk,kk)

#=========================================================================
#=================================Hexagon=================================
#=========================================================================
COEFF=2*pi/3

AA=COEFF*np.array([2,0],float)
BB=COEFF*np.array([1,sqrt(3)],float)
CC=COEFF*np.array([-1,sqrt(3)],float)
DD=COEFF*np.array([-2,0],float)
EE=COEFF*np.array([-1,-sqrt(3)],float)
FF=COEFF*np.array([1,-sqrt(3)],float)

Path=path.Path([AA,BB,CC,DD,EE,FF])
# poly = Polygon(Path)
#=========================================================================
#=================================Hexagon=================================
#=========================================================================
# COEFF=2*pi/3

kx_1=np.linspace(2,1,dim_k,endpoint=True)
ky_1=np.linspace(0,sqrt(3),dim_k,endpoint=True)

kx_2=np.linspace(1,-1,dim_k,endpoint=True)
ky_2=np.linspace(sqrt(3),sqrt(3),dim_k,endpoint=True)

kx_3=np.linspace(-1,-2,dim_k,endpoint=True)
ky_3=np.linspace(sqrt(3),0,dim_k,endpoint=True)

kx_4=np.linspace(-2,-1,dim_k,endpoint=True)
ky_4=np.linspace(0,-sqrt(3),dim_k,endpoint=True)

kx_5=np.linspace(-1,1,dim_k,endpoint=True)
ky_5=np.linspace(-sqrt(3),-sqrt(3),dim_k,endpoint=True)

kx_6=np.linspace(1,2,dim_k,endpoint=True)
ky_6=np.linspace(-sqrt(3),0,dim_k,endpoint=True)

KX=np.empty(6*dim_k,float)
KY=np.empty(6*dim_k,float)


KX[0:dim_k]=COEFF*kx_1
KX[dim_k:2*dim_k]=COEFF*kx_2
KX[2*dim_k:3*dim_k]=COEFF*kx_3
KX[3*dim_k:4*dim_k]=COEFF*kx_4
KX[4*dim_k:5*dim_k]=COEFF*kx_5
KX[5*dim_k:6*dim_k]=COEFF*kx_6

KY[0:dim_k]=COEFF*ky_1
KY[dim_k:2*dim_k]=COEFF*ky_2
KY[2*dim_k:3*dim_k]=COEFF*ky_3
KY[3*dim_k:4*dim_k]=COEFF*ky_4
KY[4*dim_k:5*dim_k]=COEFF*ky_5
KY[5*dim_k:6*dim_k]=COEFF*ky_6

#====================================================================


#====================================================================

num=0

for i in range(dim_k):
    for j in range(dim_k):
        
        Point=np.array([kx[i,j],ky[i,j]],float).reshape(1, 2)
        
        Ans=Path.contains_points(Point)
        
        if Ans==True:
            
            num+=1
            
BZ_data=np.zeros((num,2),float)

print("num = ",num)


num1=0

for i in range(dim_k):
    for j in range(dim_k):
        
        Point=np.array([kx[i,j],ky[i,j]],float).reshape(1, 2)
        
        Ans1=Path.contains_points(Point) 
        
        if Ans1==True:
            
            BZ_data[num1,0]=Point[0,0]
            BZ_data[num1,1]=Point[0,1]           
                        
            num1+=1

print("num1 = ",num1)

#=======================================================================

end_time=time.time()
Time=end_time-start_time

print("Time = ",Time)
            
#=======================================================================
 
#
np.savetxt("BZ_Data.dat",BZ_data, fmt='%.18e', delimiter=' ', newline='\n') 
#
                
print("Finished!")              
    
#=======================================================================
  
TickFontSize=35
LabelFontSize=40
TextFontSize=30

InsetFontSize=35

AxesLineWith=2
LineWidth=2

ZeroLineWidth=1

font = {'family': 'Times New Roman',
        'color':  'k',
        'weight': 'normal',
        'size': LabelFontSize
        }

rc('font',**{'family':'Times New Roman'})
rc('text', usetex=True)
plt.rcParams['text.usetex'] = True


xticks=[0,dim_k,2*dim_k,3*dim_k,4*dim_k,5*dim_k,6*dim_k]
xticklabels=[r'$\Gamma$',r'$K^\prime$',r'M',r'$\Gamma$','K','M',r'$\Gamma$']

yticks=[-2,-1,0,1,2]
yticklabels=[-2,-1,0,1,2]


Om_xticks=[0,1,2,3]
Om_xticklabels=[0,1,2,3]

Om_Nulllabels=[ ]

FSticks=[-2*pi,-pi,0,pi,2*pi]
FSticklabels=[r'$-2\pi$',r'$-\pi$','0',r'$\pi$',r'$2\pi$']
FSNullticklabels=[ ]

#============================================================================

#============================================================================


fig = plt.figure(figsize=(5,5))

plt.plot(BZ_data[:,0],BZ_data[:,1],'*',color='b')


plt.plot(AA[0],AA[1],'*',color='r')
plt.plot(BB[0],BB[1],'*',color='r')
plt.plot(CC[0],CC[1],'*',color='r')
plt.plot(DD[0],DD[1],'*',color='r')
plt.plot(EE[0],EE[1],'*',color='r')
plt.plot(FF[0],FF[1],'*',color='r')


plt.plot(KX,KY,'--',color='k',linewidth=1)


plt.xticks(FSticks,FSticklabels,fontsize=TickFontSize)
plt.yticks(FSticks,FSticklabels,fontsize=TickFontSize)

plt.xlabel(r'$k_x$', fontdict=font,fontsize=LabelFontSize)
plt.ylabel(r'$k_y$', fontdict=font,fontsize=LabelFontSize,labelpad=-4)

plt.tick_params(direction='in', length=6, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)

# plt.savefig("Band_Path.pdf",bbox_inches='tight')