# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 20:16:17 2019

@author: ajkra
"""

import numpy as np
from math import pi

import matplotlib
import matplotlib.pyplot as plt
import mplstereonet

from FUNCTIONS_coordinate_conversions import *

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('seaborn-darkgrid')
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams.update({'font.size': 16})


""" FUNCTIONS_plotting
unitspherePlot_bing_yesdata(unit_E,unit_N,unit_Z,OUT_unit_E,OUT_unit_N,OUT_unit_Z,e1,e2,e3)
unitspherePlot_bing_nodata(OUT_unit_E,OUT_unit_N,OUT_unit_Z,e1,e2,e3)
lowerspherePlot_bing(unit_E,unit_N,unit_Z,OUT_unit_E,OUT_unit_N,OUT_unit_Z,e1,e2,e3)
unitSpherePlot_vmf_yesdata(sample,unitE,unitN,unitZ,R)
unitSpherePlot_vmf_nodata(sample,R)
unitSpherePlot_vmf_yesdata_sub(ax,sample,unitE,unitN,unitZ,R)
unitSpherePlot_vmf_nodata_sub(ax,sample,R)
"""

"""
    #############################################
    ############ USE THESE ALWAYS ############### 
    #############################################

"""

def stereoplot(show,contours,title,strikes,dips,meanstrikes,meandips,numfaults,samples):
    
    plt.ioff() # Don't display figure immediately
    fig = plt.figure(figsize = (10,10))
    plt.suptitle(title)
    if numfaults == 1:
        sp = 1
        ax = fig.add_subplot(sp,sp,1,projection='stereonet')
        
        if contours == 1:
            cax = ax.density_contourf(strikes[0],dips[0],measurement='poles',cmap='Reds',sigma=5,zorder=-1)
        
        ax.pole(meanstrikes,meandips,'ro')
        ax.plane(meanstrikes,meandips,'r-',linewidth=1)
        
        for i in range(samples):      
            ax.pole(strikes[0][i],dips[0][i],'k.',markersize=1)
        
        ax.grid(True)
        
        #fig.colorbar(cax)
        
    else:
        sp = 3
        for j in range(numfaults):
            ax = fig.add_subplot(sp,sp,j+1,projection='stereonet')
            
            ax.pole(meanstrikes[j],meandips[j],'ro')
            ax.plane(meanstrikes[j],meandips[j],'r-',linewidth=1)
            
            for i in range(samples):      
                ax.pole(strikes[j][i],dips[j][i],'k.',markersize=1)
                
            ax.grid(True)
            
            faultnum = j+1
            ax.title.set_text('Fault %d \n' % faultnum)
    
    plt.legend(['Input vector (pole)','Input vector (plane)','Samples'],bbox_to_anchor=(1, 1),
           bbox_transform=plt.gcf().transFigure)
    
    ax.grid(True)
            
    if show == 1:
        plt.show()
    
    return fig


def stereoplot_eigen(show,contours,title,sample,e1,e3,mv_calc,numfaults,samples):
    
    plt.ioff() # Don't display figure immediately
    
    fig = plt.figure(figsize = (10,10))
    plt.suptitle(title)
    
    t,p,ss,a,sd = convertsamples(sample,numfaults,samples)
    e1s,e1d = convertmeans(e1,numfaults)
    mvs,mvd = convertmeans(mv_calc,numfaults)
    e3s,e3d = convertmeans(e3,numfaults)
    
    if numfaults == 1:
        subp = 1
        j = 0
        ax = fig.add_subplot(subp,subp,j+1,projection='stereonet')
        
        if contours == 1:
            ax.density_contourf(ss[j],sd[j],measurement='poles',cmap='Reds',sigma = 5,zorder=-1)
        
        ax.pole(e1s[j],e1d[j],'bo')
        ax.plane(e1s[j],e1d[j],'b-',linewidth=1)
        
        ax.pole(mvs[j],mvd[j],'ro')
        ax.plane(mvs[j],mvd[j],'r--',linewidth=1)
        
        ax.pole(e3s[j],e3d[j],'go')
        ax.plane(e3s[j],e3d[j],'g-',linewidth=1)
        
        for i in range(samples):      
            ax.pole(ss[j][i],sd[j][i],'k.',markersize=1)
            #ax.plane(vmf_strikes[j][i],vmf_dips[j][i],'k--',linewidth=0.5)
        
        ax.grid(True)
        
        
    else: 
        subp = 3
    
        for j in range(numfaults):
            ax = fig.add_subplot(subp,subp,j+1,projection='stereonet')
            
            
            ax.pole(e1s[j],e1d[j],'bo')
            ax.plane(e1s[j],e1d[j],'b-',linewidth=1)
            
            ax.pole(mvs[j],mvd[j],'ro')
            ax.plane(mvs[j],mvd[j],'r--',linewidth=1)
            
            ax.pole(e3s[j],e3d[j],'go')
            ax.plane(e3s[j],e3d[j],'g-',linewidth=1)
            
            for i in range(samples):      
                ax.pole(ss[j][i],sd[j][i],'k.',markersize=1)
                #ax.plane(vmf_strikes[j][i],vmf_dips[j][i],'k--',linewidth=0.5)
            
            ax.grid(True)
            
            faultnum = j+1
            ax.title.set_text('Fault %d \n' % faultnum)
    plt.legend(['e1 (pole)', 'e1 (plane)','Input vector (pole)', 'Input vetor (plane)', 'e3 (pole)', 'e3 (plane)',
            'Samples'],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
        
    if show == 1:
        plt.show()
    
    return fig

def stereoplot_case(show,title,sample,mv_calc,numfaults,samples):
    plt.ioff() # Don't display figure immediately
            
    fig = plt.figure(figsize = (8,8))
    plt.suptitle(title)
    
    t,p,ss,a,sd = convertsamples(sample,numfaults,samples)
    mvs,mvd = convertmeans(mv_calc,numfaults)

    if numfaults == 1:
        subp = 1
        j = 0
        ax = fig.add_subplot(subp,subp,j+1,projection='stereonet')
        
        
        ax.pole(mvs[j],mvd[j],'ro')
        ax.plane(mvs[j],mvd[j],'r--',linewidth=1)

        
        for i in range(samples):      
            ax.pole(ss[j][i],sd[j][i],'k.',markersize=1)
            #ax.plane(vmf_strikes[j][i],vmf_dips[j][i],'k--',linewidth=0.5)
        
        ax.grid(True)
        
    else: 
        subp1 = numfaults/3
        subp2 = numfaults/2
    
        for j in range(numfaults):
            ax = fig.add_subplot(subp1,subp2,j+1,projection='stereonet')
            
            ax.pole(mvs[j],mvd[j],'ro')
            ax.plane(mvs[j],mvd[j],'r--',linewidth=1)

            for i in range(samples):      
                ax.pole(ss[j][i],sd[j][i],'k.',markersize=1)
                #ax.plane(vmf_strikes[j][i],vmf_dips[j][i],'k--',linewidth=0.5)
            
            ax.grid(True)
            
            faultnum = j+1
            ax.title.set_text('Fault %d \n' % faultnum)
    plt.legend(['Fault orientation (pole)', 'Fault orientation (plane)','Samples'],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
        
    if show == 1:
        plt.show()
        
    return fig