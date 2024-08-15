#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 22:57:53 2024

@author: mehdi
"""

#=============================================================================

import time

start_time=time.time()

#==============================================================================

import numpy as np
import RPA_Susceptibility as RPAKappa

from math import pi

from parameters import *

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import rc

#==============================================================================

num_filling=9
num_phai=9

for ii in range(num_filling):
    for jj in range(num_phai):
        
        print("num_fill = ",ii)
        print("num_phai = ",jj)
    
        
        Data=np.loadtxt("Re_Bare_Kappa_pp_0_0.dat",float)
        
        dim_q=len(Data)
        print("dim_q = ",dim_q)
        
        
        print("U = ",U)
        
        
        qx=Data[:,0]
        qy=Data[:,1]
        
        #==============================================================================
        
        RPA_Kappa_00_pp=np.empty((dim_q,3),float)
        RPA_Kappa_00_pm=np.empty((dim_q,3),float)
        RPA_Kappa_00_mp=np.empty((dim_q,3),float)
        RPA_Kappa_00_mm=np.empty((dim_q,3),float)
        
        RPA_Kappa_00_intra=np.empty((dim_q,3),float)
        RPA_Kappa_00_inter=np.empty((dim_q,3),float)
        RPA_Kappa_00_total=np.empty((dim_q,3),float)
        
        
        RPA_Kappa_PM_pp=np.empty((dim_q,3),float)
        RPA_Kappa_PM_pm=np.empty((dim_q,3),float)
        RPA_Kappa_PM_mp=np.empty((dim_q,3),float)
        RPA_Kappa_PM_mm=np.empty((dim_q,3),float)
        
        RPA_Kappa_PM_intra=np.empty((dim_q,3),float)
        RPA_Kappa_PM_inter=np.empty((dim_q,3),float)
        RPA_Kappa_PM_total=np.empty((dim_q,3),float)
        
        
        RPA_Kappa_zz_pp=np.empty((dim_q,3),float)
        RPA_Kappa_zz_pm=np.empty((dim_q,3),float)
        RPA_Kappa_zz_mp=np.empty((dim_q,3),float)
        RPA_Kappa_zz_mm=np.empty((dim_q,3),float)
        
        RPA_Kappa_zz_intra=np.empty((dim_q,3),float)
        RPA_Kappa_zz_inter=np.empty((dim_q,3),float)
        RPA_Kappa_zz_total=np.empty((dim_q,3),float)
        
        
            
        
        #==============================================================================
        
        for kk in range(dim_q):
            
            #==00======================================================================
            
            RPA_Kappa_00_pp[kk,:]=RPAKappa.RPA_Kappa_00_func("Re_Bare_Kappa_pp_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_00_pm[kk,:]=RPAKappa.RPA_Kappa_00_func("Re_Bare_Kappa_pm_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_00_mp[kk,:]=RPAKappa.RPA_Kappa_00_func("Re_Bare_Kappa_mp_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_00_mm[kk,:]=RPAKappa.RPA_Kappa_00_func("Re_Bare_Kappa_mm_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            
            RPA_Kappa_00_intra[kk,0]=RPA_Kappa_00_pp[kk,0]
            RPA_Kappa_00_intra[kk,1]=RPA_Kappa_00_pp[kk,1]
            RPA_Kappa_00_intra[kk,2]=RPA_Kappa_00_pp[kk,2]+RPA_Kappa_00_mm[kk,2]
            
            RPA_Kappa_00_inter[kk,0]=RPA_Kappa_00_pp[kk,0]
            RPA_Kappa_00_inter[kk,1]=RPA_Kappa_00_pp[kk,1]
            RPA_Kappa_00_inter[kk,2]=RPA_Kappa_00_pm[kk,2]+RPA_Kappa_00_mp[kk,2]
            
            RPA_Kappa_00_total[kk,0]=RPA_Kappa_00_pp[kk,0]
            RPA_Kappa_00_total[kk,1]=RPA_Kappa_00_pp[kk,1]
            RPA_Kappa_00_total[kk,2]=RPA_Kappa_00_intra[kk,2]+RPA_Kappa_00_inter[kk,2]
            
            
            
            #==PM======================================================================
            
            RPA_Kappa_PM_pp[kk,:]=RPAKappa.RPA_Kappa_pm_func("Re_Bare_Kappa_pp_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_PM_pm[kk,:]=RPAKappa.RPA_Kappa_pm_func("Re_Bare_Kappa_pm_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_PM_mp[kk,:]=RPAKappa.RPA_Kappa_pm_func("Re_Bare_Kappa_mp_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_PM_mm[kk,:]=RPAKappa.RPA_Kappa_pm_func("Re_Bare_Kappa_mm_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            
            RPA_Kappa_PM_intra[kk,0]=RPA_Kappa_PM_pp[kk,0]
            RPA_Kappa_PM_intra[kk,1]=RPA_Kappa_PM_pp[kk,1]
            RPA_Kappa_PM_intra[kk,2]=RPA_Kappa_PM_pp[kk,2]+RPA_Kappa_PM_mm[kk,2]
            
            RPA_Kappa_PM_inter[kk,0]=RPA_Kappa_PM_pp[kk,0]
            RPA_Kappa_PM_inter[kk,1]=RPA_Kappa_PM_pp[kk,1]
            RPA_Kappa_PM_inter[kk,2]=RPA_Kappa_PM_pm[kk,2]+RPA_Kappa_PM_mp[kk,2]
            
            RPA_Kappa_PM_total[kk,0]=RPA_Kappa_PM_pp[kk,0]
            RPA_Kappa_PM_total[kk,1]=RPA_Kappa_PM_pp[kk,1]
            RPA_Kappa_PM_total[kk,2]=RPA_Kappa_PM_intra[kk,2]+RPA_Kappa_PM_inter[kk,2]
            
            
            #==zz======================================================================
             
            RPA_Kappa_zz_pp[kk,:]=RPAKappa.RPA_Kappa_zz_func("Re_Bare_Kappa_pp_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_zz_pm[kk,:]=RPAKappa.RPA_Kappa_zz_func("Re_Bare_Kappa_pm_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_zz_mp[kk,:]=RPAKappa.RPA_Kappa_zz_func("Re_Bare_Kappa_mp_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            RPA_Kappa_zz_mm[kk,:]=RPAKappa.RPA_Kappa_zz_func("Re_Bare_Kappa_mm_"+str(ii)+"_"+str(jj)+".dat",U,ii)
            
            RPA_Kappa_zz_intra[kk,0]=RPA_Kappa_zz_pp[kk,0]
            RPA_Kappa_zz_intra[kk,1]=RPA_Kappa_zz_pp[kk,1]
            RPA_Kappa_zz_intra[kk,2]=RPA_Kappa_zz_pp[kk,2]+RPA_Kappa_zz_mm[kk,2]
            
            RPA_Kappa_zz_inter[kk,0]=RPA_Kappa_zz_pp[kk,0]
            RPA_Kappa_zz_inter[kk,1]=RPA_Kappa_zz_pp[kk,1]
            RPA_Kappa_zz_inter[kk,2]=RPA_Kappa_zz_pm[kk,2]+RPA_Kappa_zz_mp[kk,2]
            
            RPA_Kappa_zz_total[kk,0]=RPA_Kappa_zz_pp[kk,0]
            RPA_Kappa_zz_total[kk,1]=RPA_Kappa_zz_pp[kk,1]
            RPA_Kappa_zz_total[kk,2]=RPA_Kappa_zz_intra[kk,2]+RPA_Kappa_zz_inter[kk,2]
            
            #==========================================================================
            
        #==============================================================================
        #
        np.savetxt("RPA_Kappa_00_total_"+str(ii)+"_"+str(jj)+".dat",RPA_Kappa_00_total, fmt='%.8e', delimiter=' ', newline='\n') 
        np.savetxt("RPA_Kappa_PM_total_"+str(ii)+"_"+str(jj)+".dat",RPA_Kappa_PM_total, fmt='%.8e', delimiter=' ', newline='\n') 
        np.savetxt("RPA_Kappa_zz_total_"+str(ii)+"_"+str(jj)+".dat",RPA_Kappa_zz_total, fmt='%.8e', delimiter=' ', newline='\n') 
        #
    
#==============================================================================

end_time=time.time()

Time=end_time-start_time

print("Time = ",Time)

#==============================================================================        


 
TicksFontSize=35
LabelFontSize=40
TextFontSize=30
TitleFontSize=32

InsetFontSize=35

AxesLineWith=2
LineWidth=2

LegendFontSize=25

ZeroLineWidth=1

font = {'family': 'Times New Roman',
        'color':  'k',
        'weight': 'normal',
        'size': LabelFontSize
        }

rc('font',**{'family':'Times New Roman'})
rc('text', usetex=True)
plt.rcParams['text.usetex'] = True


xticks=[0,dim_k_path,2*dim_k_path,3*dim_k_path-1]
xticklabels=[r'$\Gamma$',r'$K$',r'$M$',r'$\Gamma$']


yticks=[0,0.4,0.8,1.2]
yticklabels=[0,0.4,0.8,1.2]


#==========================================================================


fig=plt.figure(figsize=(10,7))  

plt.rcParams['axes.linewidth']=AxesLineWith


# ##====Raw::1========================
# plt.subplot(3, 1,1, aspect=AspectRatio)
# plt.subplots_adjust(hspace=Hspace,wspace=Wspace)

im1=plt.plot(RPA_Kappa_00_total[:,2],'g',linewidth=LineWidth,marker='^',label=r'$00$')
im2=plt.plot(RPA_Kappa_PM_total[:,2],'b',linewidth=LineWidth,label=r'$+-$')
im3=plt.plot(RPA_Kappa_zz_total[:,2],'r',linewidth=LineWidth,label=r'$zz$')

plt.axvline(x=0,ls='--',color='grey',linewidth=1)
plt.axvline(x=dim_k_path,ls='--',color='grey',linewidth=1)
plt.axvline(x=2*dim_k_path,ls='--',color='grey',linewidth=1)
plt.axvline(x=3*dim_k_path-1,ls='--',color='grey',linewidth=1)

# plt.text(30,-5,r'$r=0.6$',fontsize=TitleFontSize)  
# plt.text(x_order,y_order,r'$(a)$',color='black',fontsize=OrderFontSize)  
    

plt.xticks(xticks,xticklabels,fontsize=TicksFontSize)
plt.xlabel(r'$\bf q$', labelpad=-15,fontsize=LabelFontSize)

plt.yticks(yticks,yticklabels,fontsize=TicksFontSize)
plt.ylabel(r'$\chi^{\rm RPA}_{uu}(\bf q,0)$', labelpad=15,fontsize=LabelFontSize)

plt.legend(loc='best', bbox_to_anchor=(0.65, 0.65, 0.35, 0.35),ncols=1,fontsize=LegendFontSize)

plt.tick_params(direction='in', length=6, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)

# plt.text(-38,70,r'$\theta=0$',fontsize=TitleFontSize,rotation=90)


# plt.savefig("Fig_9.pdf",bbox_inches='tight')
# plt.savefig("Fig_9.png",bbox_inches='tight')

# ##====Raw::2========================
# plt.subplot(3, 2,3, aspect=AspectRatio)
# plt.subplots_adjust(hspace=Hspace,wspace=Wspace)

# im=plt.imshow(np.flipud(purity_fill_RD[:,:,1,2,0]),cmap=Cmap, aspect='auto')
# # cbar =plt.colorbar()
# # cbar.ax.tick_params(labelsize=ClbarFontSize)

# plt.xticks(xticks,[ ],fontsize=TicksFontSize)
# plt.xlabel(r'$ $',fontsize=LabelFontSize)

# plt.yticks(yticks,ytickslabel,fontsize=TicksFontSize)
# plt.ylabel(r'$\beta^2(\vartheta)=1-\alpha^2(\vartheta)$',fontsize=LabelFontSize)

# plt.text(-38,70,r'$\theta=\pi/4$',fontsize=TitleFontSize,rotation=90)
# plt.text(x_order,y_order,r'$(c)$',color='black',fontsize=OrderFontSize)  


# plt.subplot(3, 2,4, aspect=AspectRatio)
# plt.subplots_adjust(hspace=Hspace,wspace=Wspace)

# im=plt.imshow(np.flipud(purity_fill_RD[:,:,0,2,0]),cmap=Cmap, aspect='auto')
# # cbar =plt.colorbar()
# # cbar.ax.tick_params(labelsize=ClbarFontSize)

# plt.xticks(xticks,[ ],fontsize=TicksFontSize)
# plt.xlabel(r'$ $',fontsize=LabelFontSize)

# plt.yticks(yticks,[ ],fontsize=TicksFontSize)
# plt.ylabel(r'$ $',fontsize=LabelFontSize)

# plt.text(x_order,y_order,r'$(d)$',color='black',fontsize=OrderFontSize)  



# ##====Raw::3========================
# plt.subplot(3, 2,5, aspect=AspectRatio)
# plt.subplots_adjust(hspace=Hspace,wspace=Wspace)

# im=plt.imshow(np.flipud(purity_fill_RD[:,:,1,4,0]),cmap=Cmap, aspect='auto')
# # cbar =plt.colorbar()
# # cbar.ax.tick_params(labelsize=ClbarFontSize)

# plt.xticks(xticks,xtickslabel,fontsize=TicksFontSize)
# plt.xlabel(r'$\langle n \rangle$',fontsize=LabelFontSize)

# plt.yticks(yticks,ytickslabel,fontsize=TicksFontSize)
# plt.ylabel(r'$\beta^2(\vartheta)=1-\alpha^2(\vartheta)$',fontsize=LabelFontSize)

# plt.text(-38,70,r'$\theta=\pi/2$',fontsize=TitleFontSize,rotation=90)
# plt.text(x_order,y_order,r'$(e)$',color='black',fontsize=OrderFontSize)  


# plt.subplot(3, 2,6, aspect=AspectRatio)
# plt.subplots_adjust(hspace=Hspace,wspace=Wspace)

# im=plt.imshow(np.flipud(purity_fill_RD[:,:,0,4,0]),cmap=Cmap, aspect='auto')
# # cbar =plt.colorbar()
# # cbar.ax.tick_params(labelsize=ClbarFontSize)

# plt.xticks(xticks,xtickslabel,fontsize=TicksFontSize)
# plt.xlabel(r'$\langle n \rangle$',fontsize=LabelFontSize)

# plt.yticks(yticks,[ ],fontsize=TicksFontSize)
# plt.ylabel(r'$ $',fontsize=LabelFontSize)

# plt.text(x_order,y_order,r'$(f)$',color='black',fontsize=OrderFontSize)  


# cbar_ax = fig.add_axes([0.92,0.25, 0.02, 0.5])  # [left, bottom, width, height]
# cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
# cbar.ax.tick_params(labelsize=ClbarFontSize)
# # cbar.set_label('Colorbar Label', fontsize=LabelFontSize)


# plt.savefig("Fig_9.pdf",bbox_inches='tight')
# plt.savefig("Fig_9.png",bbox_inches='tight')




