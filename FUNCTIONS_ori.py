# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 21:25:53 2019

@author: ajkra
"""

import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from math import pi
import math
import sys

from FUNCTIONS_coordinate_conversions import *
from FUNCTIONS_plotting import *
from FUNCTIONS_base import *
from INIT_bing import *

directional = importr('Directional')

""" FUNCTIONS_bing
E,N,Z,dip,strike,unit_N,unit_E,unit_Z,length = import_bing(filename)

powertpack_def = Bing_def()
powerpack_cust = Bing_cust(lam1,lam2)

powerpack_meanvector = Bing_vec(lam1,lam2)
e1,e2,e3,eigen_matrix_r = eigenFromMean(meanvector)

orientation_matrix, orientation_matrix_list = orimat(unit_E,unit_N,unit_Z)
e1,e2,e3,eigen_values = eigenR(orientation_matrix,length)

BinghamMatrix_def = sample_bing_def(orientation_matrix,n)
BinghamMatrix_cust = sample_bing_custom(orientation_matrix,n,lam1,lam2)
BinghamMatrix_meanvector = sample_bing_vector(n,lam1,lam2,eigen_matrix_r)

BinghamMatrix_rotated = rotate_matrix(BinghamMatrix,rotation_matrix)
rotate_Hab_matrix = rotate_Hab(e1,e2,e3)
rotate_axis_matrix = rotate_axis(axis,theta)
OUT_unit_E,OUT_unit_N,OUT_unit_Z = extractUnitVectors(BinghamMatrix)
export_bing(outname,outdir,names,all_dips,all_azs)
"""

""" FUNCTION INITIALIZATION """
def import_ori(filename):
    """ Read in Filename for EJMT Tunnel mock model fault orientations """
    #filename = 'Mock_FaultOrientations.csv'
    
    """ Extract values and store in nested lists with filename identifier """
    with open(filename, 'r') as f:
        for x in range(1): # Skip header
            next(f)            
        array = np.array([[str(x) for x in line.split(',')] for line in f])
    
    names = array[:,0]
    names = names.tolist()
    
    dips = array[:,1]
    dips = dips.tolist()
    dips = [float(dips[x]) for x in range(len(dips))]
    
    azs = array[:,2]
    azs = azs.tolist()
    azs = [float(azs[x]) for x in range(len(azs))]
    
    return names,dips,azs


def orimat(unit_E,unit_N,unit_Z):
    
    orientation_matrix_list = []

    orientation_matrix_list.append(np.sum([a**2 for a in unit_E]))
    orientation_matrix_list.append(np.sum([a*b for a,b in zip(unit_E,unit_N)]))
    orientation_matrix_list.append(np.sum([a*b for a,b in zip(unit_E,unit_Z)]))
    orientation_matrix_list.append(np.sum([a*b for a,b in zip(unit_N,unit_E)]))
    orientation_matrix_list.append(np.sum([a**2 for a in unit_N]))
    orientation_matrix_list.append(np.sum([a*b for a,b in zip(unit_N,unit_Z)]))
    orientation_matrix_list.append(np.sum([a*b for a,b in zip(unit_Z,unit_E)]))
    orientation_matrix_list.append(np.sum([a*b for a,b in zip(unit_Z,unit_N)]))
    orientation_matrix_list.append(np.sum([a**2 for a in unit_Z]))
    
    orientation_matrix_vector = robjects.FloatVector(orientation_matrix_list)
    orientation_matrix = robjects.r['matrix'](orientation_matrix_vector,nrow=3)
    
    return orientation_matrix, orientation_matrix_list

eigen = robjects.r['eigen']
def eigenR(orientation_matrix, length):
    eigen_matrix = eigen(orientation_matrix)
    
    eigen_values = eigen_matrix[0]
    eigen_values = np.array(eigen_values)
    eigen_values = eigen_values/length # Normalize
    #eigen_values = [eigen_values[0]-eigen_values[2], eigen_values[1]-eigen_values[2]] # Subtract last
    
    eigen_vectors = np.array(eigen_matrix[1])
    e1 = eigen_vectors[:,0]
    e2 = eigen_vectors[:,1]
    e3 = eigen_vectors[:,2]

    return e1,e2,e3,eigen_values

def eigen_py(unitE,unitN,unitZ):
    length = len(unitE)
    matrix, matrix_list = orimat(unitE,unitN,unitZ)
    e1,e2,e3,lams = eigenR(matrix,length)
    
    return e1,e2,e3,lams



def sample_VMF(n,meanvector,kappa):
    vmf = np.array(directional.rvmf(n, robjects.FloatVector(meanvector),kappa)) # kappa = 2.022  
    return vmf

def sample_bing_def(orientation_matrix,n):
    print('\n Sampling orientation realizations... \n')
    powerpack = Bing_def()
    
    BinghamMatrix_r = powerpack.rbingham(n,orientation_matrix)
    BinghamMatrix = np.array(BinghamMatrix_r)
    
    return BinghamMatrix

def sample_bing_custom(orientation_matrix,n,lam1,lam2,lam3):
    print('\n Sampling orientation realizations... \n')
    powerpack1 = Bing_cust(lam1,lam2,lam3)
    
    BinghamMatrix_r_c = powerpack1.rbingham(n,orientation_matrix)
    BinghamMatrix_c = np.array(BinghamMatrix_r_c)
    
    return BinghamMatrix_c

def sample_kent(n,kappa,mv,beta):
    KentMatrix_r = directional.rkent(n,kappa,robjects.FloatVector(mv),beta)
    KentMatrix = np.array(KentMatrix_r)
    
    return KentMatrix
   
""" Following function rotate_axis courtesy of StackOverflow user unutbu, 2011 
    Works as intended (rotates around given axis)"""
def rotate_axis(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    A = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    
    return A

""" Apply rotation matrix to input matrix """
def rotate_matrix(BinghamMatrix,rotation_matrix):
    
    n = np.size(BinghamMatrix,0)
    
    BinghamMatrix_rotated = np.zeros((n,3))

    for i in range(len(BinghamMatrix)):

        
        """ CONSIDER USING A ROTATION ABOUT NEW E1 AXIS """
    
        BinghamMatrix_rotated[i,:] = np.transpose(np.dot(rotation_matrix,np.resize(BinghamMatrix[i],(3,1))))
        
        holder0 = BinghamMatrix_rotated[i,0]/(np.sqrt(sum([a**2 for a in BinghamMatrix_rotated[i,:]])))
        holder1 = BinghamMatrix_rotated[i,1]/(np.sqrt(sum([a**2 for a in BinghamMatrix_rotated[i,:]])))
        holder2 = BinghamMatrix_rotated[i,2]/(np.sqrt(sum([a**2 for a in BinghamMatrix_rotated[i,:]])))
        
        BinghamMatrix_rotated[i,0] = holder0
        BinghamMatrix_rotated[i,1] = holder1
        BinghamMatrix_rotated[i,2] = holder2
        
    return BinghamMatrix_rotated


""" Following function update_progress is courtesy of StackOverflow user Brian
    Khuu, https://stackoverflow.com/questions/3160699/python-progress-bar """
# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 25 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def rotate_e3_matrix(sample,e3,az,dip,meanvector):
    strike_e3 = dipazToStrike(az+90)
    trend_e3 = strikeToTrend(strike_e3)
    
    plunge = dipToPlunge(dip)
    
    e3_thresh = 0.01
    rotate_increment = 0.25*pi/180

    rotation_matrix = rotate_axis(meanvector,rotate_increment)
    rotated_sample_mat = rotate_matrix(sample,rotation_matrix)
    rotated_sample = extractUnitVectors(rotated_sample_mat)
    e1,e2,e3,lams = eigen_py(rotated_sample[0],rotated_sample[1],rotated_sample[2])
    
    for i in range(721):
        if abs(e3[2]) > e3_thresh:
            
            rotation_matrix = rotate_axis(meanvector,rotate_increment)
            rotated_sample_mat = rotate_matrix(rotated_sample_mat,rotation_matrix)
            rotated_sample = extractUnitVectors(rotated_sample_mat)
            
            e1,e2,e3,lams = eigen_py(rotated_sample[0],rotated_sample[1],rotated_sample[2])
            
        update_progress(i/720)
    
    print('Rotation aligned...')
    
    return rotated_sample_mat


    

def extractUnitVectors(BinghamMatrix):
    OUT_unit_E = BinghamMatrix[:,0]
    OUT_unit_N = BinghamMatrix[:,1]
    OUT_unit_Z = BinghamMatrix[:,2]
    
    return OUT_unit_E,OUT_unit_N,OUT_unit_Z

def extract_eigen(bing_samples,length):
    
    orientationmatrix,orimat_list = orimat(bing_samples[:,0],bing_samples[:,1],bing_samples[:,2])
    e1,e2,e3,lams = eigenR(orientationmatrix,length)
    
    return e1,e2,e3,lams
    

def export_bing(outname,names,all_dips,all_azs,numfaults,n):  
    
    with open(outname, 'w') as f:
        f.write('Fault Name,Realization Number,Dip,Dip Az\n')
        
        for i in range(numfaults):
            for j in range(n):
                f.write('%s,%f,%f,%f\n' % (names[i],j+1,all_dips[i][j],all_azs[i][j]))
    
    f.close()
    
    return
    
    

    
    
