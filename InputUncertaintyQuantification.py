# -*- coding: utf-8 -*-

"""
    ##############################################
    ############### Description ################## 
    ##############################################

From publication manuscript: 
    "Uncertainty assessment for 3D geologic modeling 
    of fault zones based on geologic inputs and prior knowledge"
    
Section 4: Input Uncertainty Characterization

      Author:   Ashton Krajnovich
    University: Colorado School of Mines (Golden, CO, USA)
    Department: Geology and Geological Engineering
    Consortium: University Transportation Center for Underground 
                Transportation Infrastructure (UTC-UTI)
                
-----------------------
| Four primary steps: |
-----------------------
                
1) Input data loading and preprocessing

2) Probability distribution characterization

3) Monte Carlo simulation

4) Data visualization and exporting

--------------------
| Four data types: |
--------------------

1) Surface trace (3D .csv piecewise-linear polyline)
    
2) Structural orientaiton (Orientation vector defined by dip angle 
    and dip azimuth)

3) Vertical termination depth (Predefined set of fixed elevations termination 
    surfaces, fault aspect ratio)
    
4) Fault zone thickness (Scalar value)
    
"""

""" 
    ##############################################
    ########## Libraries and Functions ########### 
    ##############################################
"""

# Python Standards
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from math import pi
import os

# Monte Carlo simulation
import pymc3 as pm

# Visualization
plt.style.use('seaborn-darkgrid')
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams.update({'font.size': 16})

# Bespoke functions
from FUNCTIONS_base import *


""" 
    ##############################################
    ########### General Initialization ########### 
    ##############################################
"""
import time
start = time.time()

""" Sample number """
samples = 1000
n = samples

""" Number of and names for faults to perturb """
numfaults = 1
names = ['F2']

""" Geologic model extents """
modelExtent_x = (420988,421624) 
modelExtent_y = (4391323,4392207) 
modelExtent_z = (3820,2800)

""" Directory specification """
exp_dir = "C:\\Users\\ajkra\\Documents\\UTC-UTI\\EXPERIMENTS"
exp_name = 'BenchmarkTiming'

outfolder = '%s\\%s_%d' % (exp_dir,exp_name,samples)
figoutfolder = "C:\\Users\\ajkra\\Google Drive (akrajnov@mymail.mines.edu)\\UTC-UTI\\EJMT\\Scripting\\PythonScripting\\PUBLISHING\\Geologic-Model-Input-Uncertainty-Characterization\\EX_PublishingFigures\\UpdatedFigures"#'%s\\Figures' % (outfolder)

outnames_polys = [] # automatically uses inputPolyName_samples.csv
outname_ori = '%s\\FaultOrientation_Realizations.csv' % (outfolder)
outname_term = '%s\\FaultTermination_Realizations.csv' % (outfolder)
outname_thick = '%s\\FaultThickness_Realizations.csv' % (outfolder)
poly_figname = '%s\\PolylinePerturbation' % (figoutfolder)
ori_figname = '%s\\OrientationPerturbation' % (figoutfolder)
term_figname = '%s\\TerminationPerturbation' % (figoutfolder)
thick_figname = '%s\\ThicknessPerturbation' % (figoutfolder)

""" Export and figure options """
EXPORT = 0
EXPORT_FIG = 0
figshow = 1
contourshow = 0
vmfshow = 0

if EXPORT_FIG ==1:
    if not os.path.exists(figoutfolder):
        os.makedirs(figoutfolder)

if EXPORT == 1:
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

""" 
    ##############################################
    ################## INITIALIZING ############## 
    ##############################################
"""


""" --------------------
    | Thickness inputs |
    -------------------- """
""" THICK_OBS = 0 if use displacement, 1 for observed mean/std """
thick_obs = 1

""" Observed thickness """
thick_mean = 30
thick_std = 2

""" Observed displacement """
displacement_upper = 40
displacement_lower = 30


""" ------------------------
    | Thickness parameters |
    ------------------------ """

""" Unobserved thickness, estimate from displacement with power-law 
    via Torabi et al., 2019 (y = 0.02x^0.72 )"""
    
logb_mean = np.log10(0.02)
logb_range = 0.195 #range obtained from Matlab curve fitting toolbox on digitized data
m_mean = 0.72
m_range = 0.1108 #range obtained from Matlab curve fitting toolbox on digitized data

logb_std = logb_range/3.92
m_std = m_range/3.92
    
n_crystalline = 142

CoreVsZone = 100


""" -------------------
    | Polyline inputs |
    ------------------- """

""" Read in Fault polyline filenames for EJMT Tunnel """
poly_filenames = ['F2.csv']

""" Coloring parameters """
colors = plt.cm.rainbow(np.linspace(0.3,0.85,numfaults))
colors_samples = plt.cm.rainbow(np.linspace(0.2,0.85,samples))

""" Extract locations and store in nested lists with filename identifier """
locations = []
persistences = []
avg_strikes = []
quads = []
acutes = []

for i in range(numfaults):
    locations,persistences,avg_strikes,quads,acutes = preprocess_poly(i,poly_filenames[i],names,locations,persistences,avg_strikes,quads,acutes)

""" -----------------------
    | Polyline parameters |
    ----------------------- """

map_error = 40  # Maximum geographical error of known features from the geologic map
map_std = map_error/3.92
dig_std = 3.6666 # For a 1:12,000 scale map from Zhong, 1995

# Treat thickness as a 95% confidence interval on fault trace location
thickness_std = thick_mean/CoreVsZone/3.92


""" ----------------------
    | Orientation inputs |
    ---------------------- """
infile_ori = 'Mock_FaultOrientations_PolyStrike.csv'
names = ['F2']

names_ori,dips_ori,azs_ori = import_ori(infile_ori)

strikes_ori = np.zeros((numfaults,1))
trends_ori = np.zeros((numfaults,1))
plunges_ori = np.zeros((numfaults,1))

for i in range(numfaults):
    strikes_ori[i] = dipazToStrike(azs_ori[i])
    trends_ori[i] = strikeToTrend(strikes_ori[i])
    plunges_ori[i] = dipToPlunge(dips_ori[i])
    
meanvectors = []
for i in range(numfaults):
    meanvectors.append(convertmeanvector(trends_ori[i][0],plunges_ori[i][0]))

""" --------------------------
    | Orientation parameters |
    -------------------------- """

lam1 = 200 # First eigenvalue, lambda 1
kappa = lam1*2 # For VMF, useful in simulation

shape = 0.1125 # Ratio of lambda 1 to lambda 2, determines ellipticity

lam2 = lam1*shape
lam3 = 0.0001 # Very small lambda 3

lams = [lam1,lam2,lam3]
lam_normfactor = np.sum(lams)


""" ----------------------
    | Termination inputs |
    ---------------------- """
    
""" Input the predefined termination surface depths """
domain_size = 50
z_domains = np.arange(3010,3810+1,domain_size)
numdomains = len(z_domains)

model_min = (modelExtent_x[0],modelExtent_y[0])

""" Extract the mean elevation of the fault's surface trace """
fault_elevations = []
for i in range(numfaults):
    zo = np.mean(locations[i][3])
    fault_elevations.append(zo)

""" --------------------------
    | Termination parameters |
    -------------------------- """
""" Default uniform distribution parameterization """
aspect_min = 1
aspect_max = 5

""" Optional log-normal parameterization """
aspect_logn = 1 # = 1 for using log-normal aspect ratio distribution
term_mu = np.log(2.5)
term_sd = 0.35

""" Assumed deviation of observed trace length vs. true tip-to-tip length """
persistence_std = 20 


""" 
    ##############################################
    ################### SAMPLING ################# 
    ##############################################
"""


""" ----------------------
    | Thickness sampling |
    ---------------------- """
print('\n Sampling thickness realizations... \n')

thicknesses = np.zeros((samples,numfaults))
if thick_obs == 1:
    """ Thickness perturbation from observed mean/std"""
    for i in range(numfaults):
        with pm.Model() as model:
            thick_dist = pm.Normal('Thickness',mu=thick_mean,sd=thick_std)
            
            trace = pm.sample(samples,tune=2000,cores=1)
            thickness = trace['Thickness'][0:samples]
            
            thicknesses[:,i] = thickness
            
            trace_thick = trace
            
            trace_thick_varnames = ['Thickness']

elif thick_obs == 0:
    
    displacement_samples = []
    
    for i in range(numfaults):
        
        with pm.Model() as model:
            logb_dist = pm.Normal('logb_dist',mu=logb_mean, sd=logb_std)
            m_dist = pm.Normal('m_dist',mu=m_mean,sd=m_std)

            displacement_dist = pm.Uniform('Displacement',lower = displacement_lower, upper = displacement_upper)
            thickness_dist = pm.Deterministic('Thickness',10**(m_dist*np.log10(displacement_dist)+logb_dist)*CoreVsZone)
            
            trace = pm.sample(samples, tune = 2000, cores = 1)
            
            thickness = trace['Thickness'][0:samples]
            displacement_samples.append(trace['Displacement'][0:samples])
            
            thicknesses[:,i] = thickness
            
            trace_thick = trace
            
            trace_thick_varnames = ['Displacement','Thickness']


""" ---------------------
    | Polyline sampling |
    --------------------- """
print('\n Sampling polyline realizations... \n')
    
polyline_realizations = []

for i in range(numfaults):
    current_realizations = []
        
    x_realizations = np.zeros((len(locations[i][1]),samples))
    y_realizations = np.zeros((len(locations[i][1]),samples))
    z_holder = np.zeros((len(locations[i][1]),samples))
    
    currentend = len(locations[i][1])-1
    currentx = locations[i][1]
    currenty = locations[i][2]
    
    quad = quads[i]
    acute = acutes[i]
        
    if quad == 1 or quad == 3:
        x_unc = np.cos(acute*pi/180)*thickness_std
        y_unc = np.sin(acute*pi/180)*thickness_std
        
        map_x_unc = np.cos(acute*pi/180)*map_std
        map_y_unc = np.sin(acute*pi/180)*map_std
        
        dig_x_unc = np.cos(acute*pi/180)*dig_std
        dig_y_unc = np.sin(acute*pi/180)*dig_std
        
    elif quad == 2 or quad == 4:
        x_unc = np.sin(acute*pi/180)*thickness_std
        y_unc = np.cos(acute*pi/180)*thickness_std
        
        map_x_unc = np.sin(acute*pi/180)*map_std
        map_y_unc = np.cos(acute*pi/180)*map_std
        
        dig_x_unc = np.sin(acute*pi/180)*dig_std
        dig_y_unc = np.cos(acute*pi/180)*dig_std
        
    with pm.Model():
        
        """ Uncertainty based on polyline position within fault zone (thickness) """
        x1_dist = pm.Normal('x1_dist', mu=0, sd=x_unc)
        y1_dist = pm.Normal('y1_dist', mu=0, sd=y_unc)
        
        x2_dist = pm.Normal('x2_dist', mu=0, sd=x_unc)
        y2_dist = pm.Normal('y2_dist', mu=0, sd=y_unc)
        
        """ Uncertainty based on digitization error """
        x_1_d = pm.Normal('x_1_d',mu = 0, sd = dig_x_unc)
        y_1_d = pm.Normal('y_1_d',mu = 0, sd = dig_y_unc)
        
        x_2_d = pm.Normal('x_2_d',mu = 0, sd = dig_x_unc)
        y_2_d = pm.Normal('y_2_d',mu = 0, sd = dig_y_unc)
        
        """ Direction for shift of polyline endpoints based on map error """
        dir_o_1 = pm.Uniform('dir_o_1',lower = 0, upper = pi)
        dir_o_2 = pm.Uniform('dir_o_2',lower = 0, upper = pi)
        
        """ Uncertainty based on map error """
        map_o_1 = pm.Normal('map_o_1', mu = 0, sd = map_std)
        map_o_2 = pm.Normal('map_o_2', mu = 0, sd = map_std)
        
        map_u_x_1 = pm.Deterministic('map_u_x_1', map_o_1*np.sin(dir_o_1))
        map_u_y_1 = pm.Deterministic('map_u_y_1', map_o_1*np.cos(dir_o_1))
        
        map_u_x_2 = pm.Deterministic('map_u_x_2', map_o_2*np.sin(dir_o_2))
        map_u_y_2 = pm.Deterministic('map_u_y_2', map_o_2*np.cos(dir_o_2))
        
        """ Deterministic function for combined uncertainty from all sources """
        x1_perturb = pm.Deterministic('X1 location',x1_dist+x_1_d+map_u_x_1)
        x2_perturb = pm.Deterministic('X2 location',x2_dist+x_2_d+map_u_x_2)
        
        y1_perturb = pm.Deterministic('Y1 location',y1_dist+y_1_d+map_u_y_1)
        y2_perturb = pm.Deterministic('Y2 location',y2_dist+y_2_d+map_u_y_2)
        
        """ Sample the MC model """
        trace = pm.sample(samples, tune=2000, cores=1)
        
        perturb_factor_x_1 = trace['X1 location'][0:samples]
        perturb_factor_x_2 = trace['X2 location'][0:samples]
        
        perturb_factor_y_1 = trace['Y1 location'][0:samples]
        perturb_factor_y_2 = trace['Y2 location'][0:samples]
        
        trace_poly = trace
        trace_poly_varnames = ['X1 location','X2 location','Y1 location','Y2 location']
     
    """ Apply the perturbation to the polyline endpoints and propagate along trace """
    for j in range(samples):
             
        perturb_x = np.linspace(perturb_factor_x_1[j],perturb_factor_x_2[j],len(currentx))
        
        perturb_y = np.linspace(perturb_factor_y_1[j],perturb_factor_y_2[j],len(currenty))
        
        x_realizations[:,j] = currentx+perturb_x 
        y_realizations[:,j] = currenty+perturb_y
        
    polyline_realizations.append([x_realizations,y_realizations,z_holder])


""" ------------------------
    | Orientation sampling |
    ------------------------ """

bing_samples = []
bing_mvs_calc = []
e1_calc = []
e2_calc = []
e3_calc = []

for i in range(numfaults):
    """ Initialize and sample from vMF distribution to build initial 
        orientation matrix for Bingham sampling """
    vmf_samples = []
    vmf_samples_pp = []
    meanvectors_pp = []
    
    sample = sample_VMF(n,meanvectors[i],kappa)
    vmf_samples.append(sample)

    unit_E = vmf_samples[0][:,0]
    unit_N = vmf_samples[0][:,1]
    unit_Z = vmf_samples[0][:,2]
    
    """ Build orientation matrix and compute eigenvectors """
    orientation_matrix, orientation_matrix_list = orimat(unit_E,unit_N,unit_Z)
    e1,e2,e3,eigen_values = eigenR(orientation_matrix,len(unit_E))
    lams_def = eigen_values  
    
    """ Sample the Bingham matrix """
    BinghamMatrix = sample_bing_custom(orientation_matrix,n,lam1,lam2,lam3)
    
    """ Rotate to exchange the e3 and e1 vectors (required by default) """
    rotate_Hab_matrix = rotate_axis(e3+e1,pi)
    BinghamMatrix = rotate_matrix(BinghamMatrix,rotate_Hab_matrix)
    
    """ Rotate to align the e3 great circle to be straight """
    BinghamMatrix = rotate_e3_matrix(BinghamMatrix,e3,azs_ori[i],dips_ori[i],meanvectors[i])
    
    """ Extract the unit vectors and convert to lower hemisphere, then store """
    OUT_unit_E,OUT_unit_N,OUT_unit_Z = extractUnitVectors(BinghamMatrix)
    elow,nlow,zlow = convLower(OUT_unit_E,OUT_unit_N,OUT_unit_Z)
    bing_samples.append(np.transpose(np.asarray([OUT_unit_E,OUT_unit_N,OUT_unit_Z])))
    
    """ Recompute the final eigenvectors of the sampled and rotated Bingham 
        distribution """
    orientationmatrix,orimat_list = orimat(bing_samples[i][:,0],bing_samples[i][:,1],bing_samples[i][:,2])
    e1_c,e2_c,e3_c,lams_c = eigenR(orientationmatrix,samples)
    
    e1_calc.append(e1_c)
    e2_calc.append(e2_c)
    e3_calc.append(e3_c)


""" Convert the sampled vectors to structural geology formats """
bing_trends, bing_plunges, bing_strikes, bing_azs, bing_dips = convertsamples(bing_samples,numfaults, n)


""" ------------------------
    | Termination sampling |
    ------------------------ """
print('\n Sampling termination realizations... \n')
termination_realizations = []
for i in range(numfaults):
    zo = fault_elevations[i]    
    dip = dips_ori[i]
    
    """ Parameterize the MC model """
    with pm.Model():
        
        """ Define the uniform or log-normal aspect ratio distribution """
        if aspect_logn == 1:
            aspect_dist = pm.Lognormal('Aspect ratio', mu = term_mu, sd = term_sd)
        else:
            aspect_dist = pm.Uniform('Aspect ratio',lower=aspect_min,upper=aspect_max)

        """ Define the distribution for variability of persistence """
        persistence_dist = pm.Normal('Persistence', mu = 0, sd = persistence_std)
        
        """ Deterministic distribution to compute termination depth based on 
            above distributions """
        zterm_dist = pm.Deterministic('Vertical termination',zo-np.sin(dip*pi/180)*(persistences[i]+abs(persistence_dist))/aspect_dist) 

        """ Sample the MC model and extract the relevant data """
        trace = pm.sample(samples, tune=2000, cores=1)
        persistence_samples = trace['Persistence'][0:samples]
        aspect_samples = trace['Aspect ratio'][0:samples]
        zterm_samples = trace['Vertical termination'][0:samples]
        
        trace_term = trace
        trace_term_varnames = ['Persistence','Aspect ratio','Vertical termination']
    
    termination_realizations.append(zterm_samples)


""" Convert float depths to nearest domain intervals. """
termination_export = []
for i in range(numfaults):
    current_export = []
    
    for j in range(samples):
        idx = (np.abs(z_domains-termination_realizations[i][j])).argmin()
        current = z_domains[idx]
        current_export.append(current)
        
    termination_export.append(current_export) 


""" 
    ##############################################
    ################## PLOTTING ################## 
    ##############################################
"""


""" ------------------------
    | Orientation plotting |
    ------------------------ """
    
contourshow = 0
fig = stereoplot_eigen(figshow,contourshow,'Bingham samples \n R = %d, lam1 = %.2f, lam2 = %.2f' % (samples,lam1,lam2),bing_samples,e1_calc,e3_calc,meanvectors,numfaults,samples)

""" Optional choice to show the initial vMF samples used for initialization """
if vmfshow == 1:
    fig = stereoplot_eigen(figshow,contourshow,'VMF samples',vmf_samples,e1_calc,e3_calc,meanvectors,numfaults,samples)

""" Figure export """
if EXPORT_FIG == 1:
    fig.savefig(ori_figname+'.png')
    fig.savefig(ori_figname+'.pdf')
    
contourshow = 1
fig = stereoplot_eigen(figshow,contourshow,'Bingham samples \n R = %d, lam1 = %.2f, lam2 = %.2f' % (samples,lam1,lam2),bing_samples,e1_calc,e3_calc,meanvectors,numfaults,samples)

if EXPORT_FIG == 1:
    fig.savefig(ori_figname+'_CombinedAssessment.png')
    fig.savefig(ori_figname+'_CombinedAssessment.pdf') 

""" ----------------------
    | Thickness plotting |
    ---------------------- """

fig = plt.figure(figsize=(5,10),constrained_layout=True)
for i in range(numfaults):
    ax = fig.add_subplot(1,1,i+1)
    for j in range(samples):
        if j == 0:
            ax.scatter(1,thicknesses[j,i],color = colors_samples[samples-j-1],label='Realizations')
        else:
            ax.scatter(1,thicknesses[j,i],color = colors_samples[samples-j-1],alpha=0.2,label=None)
    #ax.set_title('Samples for %s, STD = %f' % (names[i],thick_std[i]),fontsize=20)
    
    if i+1>=7:
        ax.set_xlabel(None)
    if np.mod(i+1+2,3)==0:
        ax.set_ylabel('Thickness (m)')
    ax.set_ylim((0,np.max(thicknesses[:])+2))
    if numfaults == 1:
        plt.ylabel('Thickness (m)')
        #plt.xlabel('Realization #')
    if thick_obs == 0:
        plt.title('Thickness perturbation from prior knowledge \n Displacement range: %d-%d m' % (displacement_lower,displacement_upper))
    elif thick_obs == 1:
        plt.title('Thickness perturbation from measurements \n Mean: %d m, SD: %.1f m' % (thick_mean, thick_std))
    plt.legend()
    ax.set_xlim(0,2)   
    
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
plt.show()

if EXPORT_FIG == 1:
    fig.savefig(thick_figname+'.png')
    fig.savefig(thick_figname+'.pdf')

""" Custom quality assessment plots """
thick_lims = [(20,40)]
if thick_obs == 0:
    thick_lims = [(30,40),(15,40)]
plot_post(trace_thick,trace_thick_varnames,thick_lims)
fig = plt.gcf()
if EXPORT_FIG == 1:
    fig.savefig(thick_figname+'_CombinedAssessment.png')
    fig.savefig(thick_figname+'_CombinedAssessment.pdf') 


""" ---------------------
    | Polyline plotting |
    --------------------- """

fig = plt.figure(figsize=(10,10),constrained_layout=True)
for i in range(len(polyline_realizations)):
    plt.title('Polyline perturbation \n Standard deviations - Centerline definition: %d m, Digitization: %d m, Map: %d m' 
              % (np.mean(thicknesses[:,i])/CoreVsZone/3.92,dig_std,map_std))
    plt.xlabel('Easting (m)')
    plt.ylabel('Northing (m)')
    plt.xlim(modelExtent_x)
    plt.ylim(modelExtent_y)
    ax = plt.gca()
    ax.set_aspect('equal')
    
    for j in range(samples):
        if j == samples-1:
            plt.plot(polyline_realizations[i][0][:,j],polyline_realizations[i][1][:,j],color=colors_samples[samples-j-1],linewidth=1,label='Realizations')
        else:
            plt.plot(polyline_realizations[i][0][:,j],polyline_realizations[i][1][:,j],color=colors_samples[samples-j-1],linewidth=0.5,label=None, alpha=0.1)
    plt.plot(locations[i][1],locations[i][2],c='g',linewidth=1.5,linestyle='--',label=names[i])
plt.legend()
plt.show()

if EXPORT_FIG == 1:
    fig.savefig(poly_figname+'.png')
    fig.savefig(poly_figname+'.pdf') 

""" Custom quality assessment plots """
poly_lims = [(-30,30),(-30,30),(-30,30),(-30,30)]
plot_post(trace_poly,trace_poly_varnames,poly_lims)
fig = plt.gcf()
if EXPORT_FIG == 1:
    fig.savefig(poly_figname+'_CombinedAssessment.png')
    fig.savefig(poly_figname+'_CombinedAssessment.pdf') 

""" ------------------------
    | Termination plotting |
    ------------------------ """

fig = plt.figure(figsize=(6,10),constrained_layout=True)
for i in range(numfaults):
    ax1 = fig.add_subplot(1,1,i+1)  
    ax1.axhline(np.mean(locations[i][3]),color='gold',linewidth=5,label='Fault outcrop') # Average Fault Outcrop Elevation
    #ax1.axhline(3410,color='navy',linewidth=5) # Tunnel Elevation, add to legend: 'Tunnel elevation', 

    counter = 0
    for zc in termination_realizations[i]:
        counter = counter+1
        if counter == samples:
            ax1.axhline(zc,color=colors_samples[samples-counter],linewidth=1,label='Termination depth realizations')
        else: 
            ax1.axhline(zc,color=colors_samples[samples-counter],linewidth=0.5,label=None, alpha=0.2) # Fault Termination Depths    
        
    
    unique_elements, counts_elements = np.unique(termination_export[i],return_counts=True)
    for j in range(len(unique_elements)):
        if j == 0:
            ax1.scatter(1,unique_elements[j],color = 'purple',s = counts_elements[j]/samples*150, label='Discretized termination domain',zorder=3) # Termination Domains
        else: 
            ax1.scatter(1,unique_elements[j],color = 'purple',s = counts_elements[j]/samples*150, label=None,zorder=3) # Termination Domains

    ax1.set_ylim(modelExtent_z)
    ax1.invert_yaxis()
    ax1.set_xlim(0,2)   
    
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
 
    if i+1 == 1 or i+1 == 4 or i+1 == 7:
        ax1.set_ylabel('Elevation (m)')
    
    ax1.legend()
    plt.title('Vertical termination perturbation \n Aspect ratio: %.1f-%.1f' % (aspect_min,aspect_max))

plt.show()

if EXPORT_FIG == 1:
    fig.savefig(term_figname+'.png')
    fig.savefig(term_figname+'.pdf')

""" Custom quality assessment plots """
term_lims = [(-60,60),(1,5),(2800,3550)]
plot_post(trace_term,trace_term_varnames,term_lims)
fig = plt.gcf()
if EXPORT_FIG == 1:
    fig.savefig(term_figname+'_CombinedAssessment.png')
    fig.savefig(term_figname+'_CombinedAssessment.pdf') 

""" 
    ##############################################
    ################# EXPORTING ################## 
    ##############################################
"""

if EXPORT == 1:
    
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    """ THICKNESSES """
    with open(outname_thick, 'w') as f:
        f.write('Fault Name,Realization Number,Thickness\n')
        
        for i in range(numfaults):
            for j in range(samples):
                f.write('%s,%f,%f\n' % (names[i],j+1,thicknesses[j][i]))

    
    """ POLYLINES """
    for i in range(numfaults):
        for k in range(1,samples+1):
            outnames_polys.append(['%s\Polyline_' % (outfolder) + names[i] + '_%d.csv' % (k)])
    
        
    for i in range(numfaults):
        for k in range(samples):       
            with open(outnames_polys[i*samples+k][0], 'w') as f:
                f.write('Point X,Point Y,Point Z\n')
                
                for j in range(len(polyline_realizations[i][0])):
                    f.write('%f,%f,%f\n' % (polyline_realizations[i][0][j][k],polyline_realizations[i][1][j][k],polyline_realizations[i][2][j][k]))
    
    
    """ TERMINATIONS """
    with open(outname_term, 'w') as f:
        f.write('Fault Name,Realization Number,Termination Depth\n')
        
        for i in range(numfaults):
            for j in range(samples):
                f.write('%s,%f,%f\n' % (names[i],j+1,termination_export[i][j]))
                
    """ ORIENTATIONS """
    export_bing(outname_ori,names,bing_dips,bing_azs,numfaults,n)



end = time.time()

print(end-start)