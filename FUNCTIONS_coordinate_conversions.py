# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 20:30:06 2019

@author: ajkra
"""
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from math import pi
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
import os


""" FUNCTIONS_coordinate_conversions
trend = findtrend(n,e,z)
plunge = findplunge(n,e,z)

plunge = dipToPlunge(dip)
trend = strikeToTrend(strike)
strike = dipazToStrike(az)
az = strikeToDipaz(strike)

dip,strikemarker = plungeToDip(trend,plunge)
strike = trendToStrike(trend,plunge)
strike = overturnstrike(strike,strikemarker)

unit_E_lower,unit_N_lower,unit_Z_lower = convLower(unit_E,unit_N,unit_Z)
meanvector = convertmeanvector(trend,plunge)
R_norm,R_mean = findmeanvector(unitE,unitN,unitZ)
all_trends,all_plunges,all_strikes,all_azs,all_dips = convertsamples(all_samples)
"""

""" Initialize coordinate conversion functions """
def findtrend(n,e,z):
    trend = 0
    if n > 0:
        trend = np.arctan(e/n)*180/pi
    elif n < 0:
        trend = 180+np.arctan(e/n)*180/pi
    elif n == 0 and e >= 0:
        trend = 90
    elif n == 0 and e < 0:
        trend = 270

    if trend < 0:
        trend = trend + 360
    return trend 
    
def findplunge(n,e,z):
    plunge = np.arcsin(z)*180/pi
    return plunge


def plungeToDip(trend,plunge):
    dip = 90-plunge
    strikemarker = 0
    if dip > 90:
        dip = 180-dip
        strikemarker = 1
    return dip, strikemarker

def dipToPlunge(dip):
    plunge = 90-dip
    return plunge


def trendToStrike(trend,plunge):
    if trend < 270:
        strike = trend + 90
    else:
        strike = trend-270
    return strike

def strikeToTrend(strike):
    if strike > 90:
        trend = strike-90
    else:
        trend = strike+270
    return trend


def dipazToStrike(az):
    strike = az-90
    if strike < 0:
        strike = strike+360
    return strike

def strikeToDipaz(strike):
    az = strike+90
    if az > 360:
        az = az-360
    return az


def overturnstrike(strike,strikemarker):
    if strikemarker == 1:
        strike = strike + 180
        if strike > 360:
            strike = strike-360
            
    return strike



def convLower(unit_E,unit_N,unit_Z):
    if len(np.atleast_1d(unit_Z)) == 1:
        if unit_Z > 0:
            unit_Z = -unit_Z
            unit_E = -unit_E
            unit_N = -unit_N
    else:
        for i in range(len(unit_Z)):
            if unit_Z[i] > 0:
                unit_Z[i] = -unit_Z[i]
                unit_E[i] = -unit_E[i]
                unit_N[i] = -unit_N[i]
                
    return unit_E,unit_N,unit_Z
    

def convertmeanvector(trend,plunge):
    cosa = np.cos(trend*pi/180)*np.cos(plunge*pi/180) # N
    cosb = np.sin(trend*pi/180)*np.cos(plunge*pi/180) # E
    cosg = np.sin(plunge*pi/180) # Z
    
    meanvector = [cosb, cosa, cosg] # E, N, Z
    
    return meanvector


def findmeanvector(unitE,unitN,unitZ):
    rE = np.sum(unitE)
    rN = np.sum(unitN)
    rZ = np.sum(unitZ)
    
    R = [rE,rN,rZ]
    
    n = len(unitE)
    
    rE_norm = rE/n
    rN_norm = rN/n
    rZ_norm = rZ/n
    
    R_norm = [rE_norm,rN_norm,rZ_norm]
    
    holder0 = R_norm[0]/(np.sqrt(sum([a**2 for a in R_norm[:]])))
    holder1 = R_norm[1]/(np.sqrt(sum([a**2 for a in R_norm[:]])))
    holder2 = R_norm[2]/(np.sqrt(sum([a**2 for a in R_norm[:]])))
    
    rE_mean = holder0
    rN_mean = holder1
    rZ_mean = holder2
    
    R_mean = [rE_mean,rN_mean,rZ_mean]
    
    return R_norm,R_mean


def convertsamples(all_samples,numfaults, n):
    
    all_trends = []
    all_plunges = []
    all_strikes = []
    all_azs = []
    all_dips = []
    
    for i in range(numfaults):
        trends_h = []
        plunges_h = []
        strikes_h = []
        azs_h = []
        dips_h = []
        
        for j in range(n):
            trend_c = findtrend(all_samples[i][j,1],all_samples[i][j,0],all_samples[i][j,2])
            plunge_c = findplunge(all_samples[i][j,1],all_samples[i][j,0],all_samples[i][j,2])
            
            dip_c,marker = plungeToDip(trend_c,plunge_c)
            strike_c = trendToStrike(trend_c,plunge_c)
            strike_c = overturnstrike(strike_c,marker)
    
            az_c = strikeToDipaz(strike_c)
            
            
            trends_h.append(trend_c)
            plunges_h.append(plunge_c)
            strikes_h.append(strike_c)
            azs_h.append(az_c)
            dips_h.append(dip_c)
            
        all_trends.append(trends_h)
        all_plunges.append(plunges_h)
        all_strikes.append(strikes_h)
        all_azs.append(azs_h)
        all_dips.append(dips_h)

    return all_trends,all_plunges,all_strikes,all_azs,all_dips
    
def convertdata(data,numfaults, n):
    all_trends = []
    all_plunges = []
    all_strikes = []
    all_azs = []
    all_dips = []
    
    for j in range(n):
        trend_c = findtrend(data[j,1],data[j,0],data[j,2])
        plunge_c = findplunge(data[j,1],data[j,0],data[j,2])
        
        dip_c,marker = plungeToDip(trend_c,plunge_c)
        strike_c = trendToStrike(trend_c,plunge_c)
        strike_c = overturnstrike(strike_c,marker)

        az_c = strikeToDipaz(strike_c)
        
            
        all_trends.append(trend_c)
        all_plunges.append(plunge_c)
        all_strikes.append(strike_c)
        all_azs.append(az_c)
        all_dips.append(dip_c)

    return all_trends,all_plunges,all_strikes,all_azs,all_dips

def convertmeans(meanvectors,numfaults):
    
    all_strikes = []
    all_dips = []
    
    for j in range(numfaults):

        trend_c = findtrend(meanvectors[j][1],meanvectors[j][0],meanvectors[j][2])
        plunge_c = findplunge(meanvectors[j][1],meanvectors[j][0],meanvectors[j][2])
        
        dip_c,marker = plungeToDip(trend_c,plunge_c)
        strike_c = trendToStrike(trend_c,plunge_c)
        strike_c = overturnstrike(strike_c,marker)

        all_strikes.append(strike_c)
        all_dips.append(dip_c)

    return all_strikes,all_dips


def convertmean(meanvector):

    trend_c = findtrend(meanvector[1],meanvector[0],meanvector[2])
    plunge_c = findplunge(meanvector[1],meanvector[0],meanvector[2])
    
    dip_c,marker = plungeToDip(trend_c,plunge_c)
    strike_c = trendToStrike(trend_c,plunge_c)
    strike_c = overturnstrike(strike_c,marker)

    return strike_c,dip_c
    

def anglebetween(v1,v2):
    
    v1_l = v1/np.linalg.norm(v1)
    v2_l = v2/np.linalg.norm(v2)
   
    aa = np.dot(v1_l,v2_l)
    if abs(abs(aa)-1) < 0.0001:
        aa = 1
        
    angle_between = np.arccos(aa)*180/pi

    angle_between = np.min([180-angle_between,angle_between])
    
    return angle_between
    

def colvec(vec):
    col = np.reshape(np.asarray(vec),(len(vec),1))
    return col