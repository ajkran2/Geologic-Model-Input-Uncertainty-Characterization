# -*- coding: utf-8 -*-
"""
Collection of miscellaneous bespoke functions for InputUncertaintyQuantification.py
"""

import arviz
from FUNCTIONS_ori import *

def preprocess_poly(i,filename,names,locations,persistences,avg_strikes,quads,acutes):
    
    with open(filename, 'r') as f:
        for x in range(2): # Skip header
            next(f)            
        array = [[str(x) for x in line.split(',')] for line in f]   
    array = array[0:len(array)-1]
    
    x = []
    y = []
    z = []
    
    for j in range(len(array)):
        x.append(float(array[j][0]))
        y.append(float(array[j][1]))
        z.append(float(array[j][2]))   
    
    """ Compute Persistence """
    difference = np.diff([x,y,z])
    persistence = np.sum(np.sqrt(np.sum(difference*difference,0)))
    
    persistences.append(persistence)
    
    
    """ Compute Average Strike """
    xx = np.array(x)
    yy = np.array(y)
    
    yx = np.diff([xx,yy])
    segment_distance = np.sqrt(np.sum(yx*yx,0))
    segment_strikes = np.arctan2(yx[0],yx[1])*180/pi
    avg_strike = np.average(segment_strikes,axis = None, weights = segment_distance)
    if avg_strike < 0:
        avg_strike = avg_strike+360

    avg_strikes.append(avg_strike)

    locations.append([names[i],x,y,z,persistence,avg_strike])
    
    
    avg_strike = avg_strikes[i]
    
    quad = 0
    
    if 0 <= avg_strike < 90:
        quad = 1
        if avg_strike == 0:
            avg_strike = 0.01
    elif 90 <= avg_strike < 180:
        quad = 2
        if avg_strike == 90:
            avg_strike = 90.01
    elif 180 <= avg_strike < 270:
        quad = 3
        if avg_strike == 180:
            avg_strike = 180.01
    elif 270 <= avg_strike < 360:
        quad = 4
        if avg_strike == 270:
            avg_strike == 270.01
    
    acute = np.mod(avg_strike,90)
    if acute > 45:
        acute = 90-acute
    
    quads.append(quad)
    acutes.append(acute)
        
    return locations,persistences,avg_strikes,quads,acutes

def plot_trace_post(data,var_names,lims):
    n = len(var_names)
    fig,axes = plt.subplots(nrows=n,ncols=1,constrained_layout=True)
    for i in range(n):
        if n>1:
            ax = axes[i]
        else:
            ax = axes
        ax.set_xlim(lims[i])
        arviz.plot_posterior(data,var_names=var_names[i],ax=ax,bins=100,kind='hist',credible_interval=0.95)
    
    axes = arviz.plot_trace(data,var_names=var_names)
    for i in range(n):
        ax = axes[i][1]
        ax.set_ylim(lims[i])
        
    plt.show()
        