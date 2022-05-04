import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob


# change the default font family
plt.rcParams.update({'font.family':'Arial'})

def myexponential(x,A,l):
    return A*np.exp(-x/l)+1

def myexponential_shifted(x,A,l):
    return A*np.exp(-(x-2)/l)+1

def mypowerlaw(x,A,l):
    return A*x**(-l)+1

def rmsd(x,y):
    Nx = len(x)
    Ny = len(y)
    if Nx!=Ny:
        print('WARNING! Nx!=Ny. Could not calculate rmsd value')
        return 'WARNING! Nx!=Ny. Could not calculate rmsd value'
    delta = 0
    for i in range(Nx):
        delta += (x[i]-y[i])*(x[i]-y[i])
    delta = np.sqrt(delta/(Nx-1))
    return delta

Nintervals = 10
limitsmall = 11

# Input parameters for file selection:
psigma   = 1
damp     = 10
# Input booleans for file selection:
bulkdiffusion = False
substrate     = False
moresigmas    = False
big           = False
bulk_cut      = False
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
old_bulk      = False
endlocation = '/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
if moresigmas==True:
    endlocation = '/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'
endlocation_nostiff = endlocation + 'Nostiffness/'
endlocation_f = '/Diffusion_forest/D_vs_d/Nocut/'

basepath_base      = '/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = '/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'
if old_bulk==False:
    bulkfilename  = bulklocation + 'diffusion_bulk'+filestext+'_new.txt'

## Files to read
brushfilename_dyn      = endlocation + 'D_vs_d_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_f        = endlocation_f + 'D_vs_d_forest_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_stat     = endlocation_static + 'D_vs_d_static_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_nostiff  = endlocation_nostiff + 'D_vs_d_better_rms_Nestimates%i.txt' % Nintervals
## Files to write to
if big==False:
    plotname     = endlocation_static+'Dd_dyn_vs_stat_noDR_better_rms_Nestimates%i_all.png' % Nintervals
    plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut_noDR_better_rms_Nestimates%i_all.png' % Nintervals
    plotname_twoinone = endlocation_static+'Dd_dyn_vs_stat_noDR_twoinone_better_rms_Nestimates%i_all.png' % Nintervals
else:
    plotname     = endlocation_static+'Dd_dyn_vs_stat_big_noDR_better_rms_Nestimates%i_all.png' % Nintervals
    plotname_cut = endlocation_static+'Dd_dyn_vs_stat_cut_big_noDR_better_rms_Nestimates%i_all.png' % Nintervals
    plotname_twoinone = endlocation_static+'Dd_dyn_vs_stat_big_noDR_twoinone_better_rms_Nestimates%i_all.png' % Nintervals
plotname_Dz_fraction_dyn = endlocation_static+'D_vs_d_perp_div_parallel_dyn.png'
plotname_Dz_fraction_stat = endlocation_static+'D_vs_d_perp_div_parallel_stat.png'
plotname_Dz_fraction_forest = endlocation_static+'D_vs_d_perp_div_parallel_forest.png'
plotname_Dz_fraction_nostiff = endlocation_static+'D_vs_d_perp_div_parallel_nostiff.png'

# Dynamic sims:
brushfile_dyn = open(brushfilename_dyn, 'r')

lines = brushfile_dyn.readlines()
N_dyn = len(lines)-1

# ds
spacings_dyn = np.zeros(N_dyn)
# Ds
DRs_dyn = np.zeros(N_dyn)
Dxs_dyn = np.zeros(N_dyn)
Dys_dyn = np.zeros(N_dyn)
Dzs_dyn = np.zeros(N_dyn)
Dparallel_dyn = np.zeros(N_dyn)
# Ds, stdv
DRs_stdv_dyn = np.zeros(N_dyn)
Dxs_stdv_dyn = np.zeros(N_dyn)
Dys_stdv_dyn = np.zeros(N_dyn)
Dzs_stdv_dyn = np.zeros(N_dyn)
Dparallel_stdv_dyn = np.zeros(N_dyn)

for i in range(1,N_dyn+1):
    words = lines[i].split()
    j = i-1
    
    spacings_dyn[j] = float(words[0])
    DRs_dyn[j] = float(words[1])
    Dzs_dyn[j] = float(words[3])
    Dparallel_dyn[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_dyn[j] = float(words[2])
    Dzs_stdv_dyn[j] = float(words[4])
    Dparallel_stdv_dyn[j] = float(words[6])
    
brushfile_dyn.close()

##########

# Static sims:
brushfile_stat = open(brushfilename_stat, 'r')

lines = brushfile_stat.readlines()
N_stat = len(lines)-1

# ds
spacings_stat = np.zeros(N_stat)
# Ds
DRs_stat = np.zeros(N_stat)
Dxs_stat = np.zeros(N_stat)
Dys_stat = np.zeros(N_stat)
Dzs_stat = np.zeros(N_stat)
Dparallel_stat = np.zeros(N_stat)
# Ds, stdv
DRs_stdv_stat = np.zeros(N_stat)
Dxs_stdv_stat = np.zeros(N_stat)
Dys_stdv_stat = np.zeros(N_stat)
Dzs_stdv_stat = np.zeros(N_stat)
Dparallel_stdv_stat = np.zeros(N_stat)

for i in range(1,N_stat+1):
    words = lines[i].split()
    j = i-1
    
    spacings_stat[j] = float(words[0])
    DRs_stat[j] = float(words[1])
    Dzs_stat[j] = float(words[3])
    Dparallel_stat[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_stat[j] = float(words[2])
    Dzs_stdv_stat[j] = float(words[4])
    Dparallel_stdv_stat[j] = float(words[6])
    
brushfile_stat.close()

#######

# Forest sims:
brushfile_f = open(brushfilename_f, 'r')

lines = brushfile_f.readlines()
N_f = len(lines)-1

# ds
spacings_f = np.zeros(N_f)
# Ds
DRs_f = np.zeros(N_f)
Dxs_f = np.zeros(N_f)
Dys_f = np.zeros(N_f)
Dzs_f = np.zeros(N_f)
Dparallel_f = np.zeros(N_f)
# Ds, stdv
DRs_stdv_f = np.zeros(N_f)
Dxs_stdv_f = np.zeros(N_f)
Dys_stdv_f = np.zeros(N_f)
Dzs_stdv_f = np.zeros(N_f)
Dparallel_stdv_f = np.zeros(N_f)

for i in range(1,N_f+1):
    words = lines[i].split()
    j = i-1
    
    spacings_f[j] = float(words[0])
    DRs_f[j] = float(words[1])
    Dzs_f[j] = float(words[3])
    Dparallel_f[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_f[j] = float(words[2])
    Dzs_stdv_f[j] = float(words[4])
    Dparallel_stdv_f[j] = float(words[6])
    
brushfile_f.close()

#####
# Nostiff sims:
brushfile_nostiff = open(brushfilename_nostiff, 'r')

lines = brushfile_nostiff.readlines()
N_nostiff = len(lines)-1

# ds
spacings_nostiff = np.zeros(N_nostiff)
# Ds
DRs_nostiff = np.zeros(N_nostiff)
Dxs_nostiff = np.zeros(N_nostiff)
Dys_nostiff = np.zeros(N_nostiff)
Dzs_nostiff = np.zeros(N_nostiff)
Dparallel_nostiff = np.zeros(N_nostiff)
# Ds, stdv
DRs_stdv_nostiff = np.zeros(N_nostiff)
Dxs_stdv_nostiff = np.zeros(N_nostiff)
Dys_stdv_nostiff = np.zeros(N_nostiff)
Dzs_stdv_nostiff = np.zeros(N_nostiff)
Dparallel_stdv_nostiff = np.zeros(N_nostiff)

for i in range(1,N_nostiff+1):
    words = lines[i].split()
    j = i-1
    
    spacings_nostiff[j] = float(words[0])
    DRs_nostiff[j] = float(words[1])
    Dzs_nostiff[j] = float(words[3])
    Dparallel_nostiff[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_nostiff[j] = float(words[2])
    Dzs_stdv_nostiff[j] = float(words[4])
    Dparallel_stdv_nostiff[j] = float(words[6])
    
brushfile_nostiff.close()


###
#Bulk:

bulkfile  = open(bulkfilename, 'r')
# D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
if old_bulk==True: # Worse rms
    bulklines = bulkfile.readlines()
    bulkline  = bulklines[1]
    words     = bulkline.split()
    
    # Ds
    DRs_bulk = float(words[0])
    Dzs_bulk = float(words[4])
    Dparallel_bulk = float(words[8])
    
    # Ds, stdv
    DRs_stdv_bulk = float(words[1])
    Dzs_stdv_bulk = float(words[5])
    Dparallel_stdv_bulk = float(words[9])
else:
    bulklines = bulkfile.readlines()
    bulkline  = bulklines[1]
    words     = bulkline.split()
    
    # Ds
    DRs_bulk = float(words[1])
    Dzs_bulk = float(words[3])
    Dparallel_bulk = float(words[5])
    
    # Ds, stdv
    DRs_stdv_bulk = float(words[2])
    Dzs_stdv_bulk = float(words[4])
    Dparallel_stdv_bulk = float(words[6])
    
bulkfile.close()

# Divide by bulk:
for i in range(N_stat):
    DRnew = DRs_stat[i]/DRs_bulk
    Dznew = Dzs_stat[i]/DRs_bulk
    Dparnew = Dparallel_stat[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_stat[i] = abs(DRnew)*np.sqrt((DRs_stdv_stat[i]/DRs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_stat[i] = abs(Dznew)*np.sqrt((Dzs_stdv_stat[i]/Dzs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_stat[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_stat[i]/Dparallel_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_stat[i] = DRnew
    Dzs_stat[i] = Dznew
    Dparallel_stat[i] = Dparnew

for i in range(N_dyn):
    DRnew = DRs_dyn[i]/DRs_bulk
    Dznew = Dzs_dyn[i]/DRs_bulk
    Dparnew = Dparallel_dyn[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_dyn[i] = abs(DRnew)*np.sqrt((DRs_stdv_dyn[i]/DRs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_dyn[i] = abs(Dznew)*np.sqrt((Dzs_stdv_dyn[i]/Dzs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_dyn[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_dyn[i]/Dparallel_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_dyn[i] = DRnew
    Dzs_dyn[i] = Dznew
    Dparallel_dyn[i] =  Dparnew

####
for i in range(N_f):
    DRnew = DRs_f[i]/DRs_bulk
    Dznew = Dzs_f[i]/DRs_bulk
    Dparnew = Dparallel_f[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_f[i] = abs(DRnew)*np.sqrt((DRs_stdv_f[i]/DRs_f[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_f[i] = abs(Dznew)*np.sqrt((Dzs_stdv_f[i]/Dzs_f[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_f[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_f[i]/Dparallel_f[i])**2 +(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_f[i] = DRnew
    Dzs_f[i] = Dznew
    Dparallel_f[i] = Dparnew

for i in range(N_nostiff):
    DRnew = DRs_nostiff[i]/DRs_bulk
    Dznew = Dzs_nostiff[i]/DRs_bulk
    Dparnew = Dparallel_nostiff[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_nostiff[i] = abs(DRnew)*np.sqrt((DRs_stdv_nostiff[i]/DRs_nostiff[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_nostiff[i] = abs(Dznew)*np.sqrt((Dzs_stdv_nostiff[i]/Dzs_nostiff[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_nostiff[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_nostiff[i]/Dparallel_nostiff[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_nostiff[i] = DRnew
    Dzs_nostiff[i] = Dznew
    Dparallel_nostiff[i] = Dparnew

minforsmall = 0
maxforsmall = 0
i = 0
d = 1
while d<limitsmall:
    ## Max vals:
    Dthismax_zdyn = Dzs_dyn[i]+Dzs_stdv_dyn[i]
    Dthismax_pdyn = Dparallel_dyn[i]+Dparallel_stdv_dyn[i]
    Dthismax_zstat = Dzs_stat[i]+Dzs_stdv_stat[i]
    Dthismax_pstat = Dparallel_stat[i]+Dparallel_stdv_stat[i]
    ## Min vals:
    Dthismin_zdyn = Dzs_dyn[i]-Dzs_stdv_dyn[i]
    Dthismin_pdyn = Dparallel_dyn[i]-Dparallel_stdv_dyn[i]
    Dthismin_zstat = Dzs_stat[i]-Dzs_stdv_stat[i]
    Dthismin_pstat = Dparallel_stat[i]-Dparallel_stdv_stat[i]
    ###### Testing if larger:
    if Dthismax_zdyn>maxforsmall:
        maxforsmall=Dthismax_zdyn
    if Dthismax_pdyn>maxforsmall:
        maxforsmall=Dthismax_pdyn
    if Dthismax_zstat>maxforsmall:
        maxforsmall=Dthismax_zstat
    if Dthismax_pstat>maxforsmall:
        maxforsmall=Dthismax_pstat
    ##### Testing if smaller:
    if Dthismin_zdyn<minforsmall:
        minforsmall=Dthismin_zdyn
    if Dthismin_pdyn<minforsmall:
        minforsmall=Dthismin_pdyn
    if Dthismin_zstat<minforsmall:
        minforsmall=Dthismin_zstat
    if Dthismin_pstat<minforsmall:
        minforsmall=Dthismin_pstat
    i+=1
    d = spacings_dyn[i] 

# Double check for forest:
i = 0
d = 1
while d<limitsmall:
    # Max vals:
    Dthismax_zf = Dzs_f[i]+Dzs_stdv_f[i]
    Dthismax_pf = Dparallel_f[i]+Dparallel_stdv_f[i]
    # Min vals:
    Dthismin_zf = Dzs_f[i]-Dzs_stdv_f[i]
    Dthismin_pf = Dparallel_f[i]-Dparallel_stdv_f[i]
    #### Testing if larger:
    if Dthismax_zf>maxforsmall:
        maxforsmall=Dthismax_zf
    if Dthismax_pf>maxforsmall:
        maxforsmall=Dthismax_pf
    #### Testing if smaller:
    if Dthismin_zf<minforsmall:
        minforsmall=Dthismin_zf
    if Dthismin_pf<minforsmall:
        minforsmall=Dthismin_pf
    i+=1
    d = spacings_f[i] 


# Double check for nostiff
i = 0
d = 1
while d<11:
    # Max vals:
    Dthismax_zn = Dzs_nostiff[i]+Dzs_stdv_nostiff[i]
    Dthismax_pn = Dparallel_nostiff[i]+Dparallel_stdv_nostiff[i]
    # Min vals:
    Dthismin_zn = Dzs_nostiff[i]-Dzs_stdv_nostiff[i]
    Dthismin_pn = Dparallel_nostiff[i]-Dparallel_stdv_nostiff[i]
    #### Testing if larger:
    if Dthismax_zn>maxforsmall:
        maxforsmall=Dthismax_zn
    if Dthismax_pn>maxforsmall:
        maxforsmall=Dthismax_pn
    #### Testing if smaller:
    if Dthismin_zn<minforsmall:
        minforsmall=Dthismin_zn
    if Dthismin_pn<minforsmall:
        minforsmall=Dthismin_pn
    i+=1
    d = spacings_nostiff[i] 



if minforsmall<0:
    minforsmall*=1.2
else:
    minforsmall*=0.8
maxforsmall*=1.05
print('maxforsmall:',maxforsmall)

'''
if big==False:
    plt.figure(figsize=(8,5))
    ax = plt.subplot(111)    
    # Data points
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.plot(spacings_f, Dzs_f, '-v' , color='r',label=r'$D_\perp$, straight')
    ax.plot(spacings_f, Dparallel_f, '-1', color='black', label=r'$D_\parallel$, straight')
    ax.plot(spacings_nostiff, Dzs_nostiff, '-X', color='lightcoral', label=r'$D_\perp$, no stiffness')
    ax.plot(spacings_nostiff, Dparallel_nostiff, '-<', color='lightgray',label=r'$D_\parallel$, no stiffness')
    # Fill between
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    ax.fill_between(spacings_f, Dzs_f+Dzs_stdv_f, Dzs_f-Dzs_stdv_f, color='r',alpha=0.2)
    ax.fill_between(spacings_f, Dparallel_f+Dparallel_stdv_f, Dparallel_f-Dparallel_stdv_f, color='black',alpha=0.2)
    ax.fill_between(spacings_nostiff, Dzs_nostiff+Dzs_stdv_nostiff, Dzs_nostiff-Dzs_stdv_nostiff,color='lightcoral',alpha=0.2)
    ax.fill_between(spacings_nostiff, Dparallel_nostiff+Dparallel_stdv_nostiff, Dparallel_nostiff-Dparallel_stdv_nostiff, color='lightgray',alpha=0.5)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$')
    else:
        plt.xlabel(r'$d$ (nm)')
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
else:
    plt.figure(figsize=(16,10))
    plt.rc('xtick', labelsize=14) 
    plt.rc('ytick', labelsize=14) 
    ax = plt.subplot(111)
    # Data points
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.plot(spacings_f, Dzs_f, '-v' , color='r',label=r'$D_\perp$, straight')
    ax.plot(spacings_f, Dparallel_f, '-1', color='black', label=r'$D_\parallel$, straight')
    ax.plot(spacings_nostiff, Dzs_nostiff, '-X', color='lightcoral', label=r'$D_\perp$, no stiffness')
    ax.plot(spacings_nostiff, Dparallel_nostiff, '-<', color='lightgray',label=r'$D_\parallel$, no stiffness')
    # Fill between:
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$', fontsize=14)
    else:
        plt.xlabel(r'$d$ (nm)', fontsize=14)
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., prop={'size': 20})
plt.savefig(plotname)

if big==False:
    plt.figure(figsize=(6.4,5))
    ax = plt.subplot(111)    
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$')
    else:
        plt.xlabel(r'$d$ (nm)')
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
else:
    plt.figure(figsize=(12.8,10))
    ax = plt.subplot(111)
    ax.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    if moresigmas==True:
        plt.xlabel(r'$d/\sigma_b$', fontsize=14)
    else:
        plt.xlabel(r'$d$ (nm)', fontsize=14)
    plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
    ax.axis([0,11,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(plotname_cut)

plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
################## twoinone ####################################
Nattempt = len(spacings_dyn) # Or 100, maybe? Need to add another line, in that case.
attempt = np.zeros(Nattempt)
for i in range(Nattempt):
    attempt[i] = 1-np.pi*(1/spacings_dyn[i])**2
if big==False:
    print("Plotting it")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4),dpi=300)
    ax1.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax1.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax1.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax1.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    ax1.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax1.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax1.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax1.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        ax1.set_xlabel(r'$d/\sigma_b$', fontsize=13)
        ax1.set_ylabel('$D/D_{\mathregular{bulk}}$', fontsize=13)
    else:
        ax1.set_xlabel(r'$d$ (nm)', fontsize=13)
        ax1.set_ylabel('$D/D_{\mathregular{bulk}}$', fontsize=13)
    ax1.set_title('A',loc='left',fontsize=14)
    ax1.legend(loc="lower right")#,fontsize=10)
    # Fill between
    ax2.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax2.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax2.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax2.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    ax2.fill_between(spacings_f, Dzs_f+Dzs_stdv_f, Dzs_f-Dzs_stdv_f, color='r',alpha=0.2)
    ax2.fill_between(spacings_f, Dparallel_f+Dparallel_stdv_f, Dparallel_f-Dparallel_stdv_f, color='black',alpha=0.2)
    ax2.fill_between(spacings_nostiff, Dzs_nostiff+Dzs_stdv_nostiff, Dzs_nostiff-Dzs_stdv_nostiff,color='lightcoral',alpha=0.2)
    ax2.fill_between(spacings_nostiff, Dparallel_nostiff+Dparallel_stdv_nostiff, Dparallel_nostiff-Dparallel_stdv_nostiff, color='lightgray',alpha=0.5)
    # Plot data
    #ax2.plot(spacings_dyn, attempt, '-o', label=r'$1-\pi d^{-2}$')
    ax2.plot(spacings_dyn, Dzs_dyn, '-o' , color='b')#,label=r'$D_\perp$, dyn.')
    ax2.plot(spacings_dyn, Dparallel_dyn, '-*', color='g')#, label=r'$D_\parallel$, dyn.')
    ax2.plot(spacings_stat, Dzs_stat, '-s', color='c')#, label=r'$D_\perp$, stat.')
    ax2.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen')#,label=r'$D_\parallel$, stat.')
    ax2.plot(spacings_f, Dzs_f, '-v' , color='r',label=r'$D_\perp$, straight')
    ax2.plot(spacings_f, Dparallel_f, '-1', color='black', label=r'$D_\parallel$, straight')
    ax2.plot(spacings_nostiff, Dzs_nostiff, '-X', color='lightcoral', label=r'$D_\perp$, no stiff.')
    ax2.plot(spacings_nostiff, Dparallel_nostiff, '-<', color='lightgray',label=r'$D_\parallel$, no stiff.')
    if moresigmas==True:
        ax2.set_xlabel(r'$d/\sigma_b$', fontsize=13)
        ax2.set_ylabel(ylabel=r'$D/D_{\mathregular{bulk}}$', fontsize=13)
    else:
        ax2.set_xlabel(r'$d$ (nm)', fontsize=13)
        ax2.set_ylabel(ylabel=r'$D/D_{\mathregular{bulk}}$', fontsize=13)
    ax2.axis([0,limitsmall,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    ax2.set_title('B',loc='left',fontsize=14)
    ax2.legend(loc="upper left")#,fontsize=10)
    fig.tight_layout()
    plt.show()
    fig.savefig(plotname_twoinone)
else:
    print("Plottin it")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4),dpi=300)
    ax1.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax1.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax1.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax1.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    #ax1.plot(spacings_dyn, attempt, '-o' ,label=r'$1-\pi d^{-2}$')
    ax1.plot(spacings_dyn, Dzs_dyn, '-o' , color='b',label=r'$D_\perp$, dyn.')
    ax1.plot(spacings_dyn, Dparallel_dyn, '-*', color='g', label=r'$D_\parallel$, dyn.')
    ax1.plot(spacings_stat, Dzs_stat, '-s', color='c', label=r'$D_\perp$, stat.')
    ax1.plot(spacings_stat, Dparallel_stat, '-d', color='limegreen',label=r'$D_\parallel$, stat.')
    if moresigmas==True:
        ax1.set_xlabel(r'$d/\sigma_b$', fontsize=13)
        ax1.set_ylabel('$D/D_{\mathregular{bulk}}$', fontsize=13)
    else:
        ax1.set_xlabel(r'$d$ (nm)', fontsize=13)
        ax1.set_ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=13)
    ax1.set_title('A',loc='left',fontsize=14)
    # Fill between
    ax2.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='b',alpha=0.2)
    ax2.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
    ax2.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='c', alpha=0.2)
    ax2.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen',alpha=0.2)
    ax2.fill_between(spacings_f, Dzs_f+Dzs_stdv_f, Dzs_f-Dzs_stdv_f, color='r',alpha=0.2)
    ax2.fill_between(spacings_f, Dparallel_f+Dparallel_stdv_f, Dparallel_f-Dparallel_stdv_f, color='black',alpha=0.2)
    ax2.fill_between(spacings_nostiff, Dzs_nostiff+Dzs_stdv_nostiff, Dzs_nostiff-Dzs_stdv_nostiff,color='lightcoral',alpha=0.2)
    ax2.fill_between(spacings_nostiff, Dparallel_nostiff+Dparallel_stdv_nostiff, Dparallel_nostiff-Dparallel_stdv_nostiff, color='lightgray',alpha=0.5)
    # Plot data
    #ax2.plot(spacings_dyn, attempt, '-o' ,label=r'$1-\pi d^{-2}$')
    ax2.plot(spacings_dyn, Dzs_dyn, '-o' , color='b')
    ax2.plot(spacings_dyn, Dparallel_dyn, '-*', color='g')
    ax2.plot(spacings_stat, Dzs_stat, '-s', color='c')
    ax2.plot(spacings_stat, Dparallel_stat, '-d')#, color='limegreen',label=r'$D_\parallel$, stat.')
    ax2.plot(spacings_f, Dzs_f, '-v' , color='r',label=r'$D_\perp$, straight')
    ax2.plot(spacings_f, Dparallel_f, '-1', color='black', label=r'$D_\parallel$, straight')
    ax2.plot(spacings_nostiff, Dzs_nostiff, '-X', color='lightcoral', label=r'$D_\perp$, no stiff.')
    ax2.plot(spacings_nostiff, Dparallel_nostiff, '-<', color='lightgray',label=r'$D_\parallel$, no stiff.')
    if moresigmas==True:
        ax2.set_xlabel(r'$d/\sigma_b$', fontsize=13)
        ax2.set_ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=13)
    else:
        ax2.set_xlabel(r'$d$ (nm)', fontsize=13)
        ax2.set_ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=13)
    ax2.axis([0,limitsmall,minforsmall,maxforsmall]) # 6e-7 before we divided by Dbulk
    #ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.set_title('B',loc='left',fontsize=14)
    ax2.legend(loc="upper left")
    fig.tight_layout()
    plt.show()
    fig.savefig(plotname_twoinone)

### Dyn, perp vs parallel

spacings_fraction = []
Dzs_fraction      = []
Dzs_diff          = []
Dzs_rms_fraction  = []
Dzs_rms_diff      = []
ones              = []

spacings_few = []
Dzs_fraction_few = []
Dzs_fraction_rms_few = []

spacings_plot = []
Dzs_fraction_plot = []
Dzs_fraction_rms_plot = []

i = 0
spacing = spacings_f[i]
while spacing<=25:
    spacings_fraction.append(spacing)    
    Dznew = Dzs_dyn[i]/Dparallel_dyn[i]
    # Ds, stdv
    Dzs_stdv = abs(Dznew)*np.sqrt((Dzs_stdv_dyn[i]/Dzs_dyn[i])**2+(Dparallel_stdv_dyn[i]/Dparallel_dyn[i])**2)
    Dzs_stdv_diff = np.sqrt((Dzs_stdv_dyn[i])**2+(Dparallel_stdv_dyn[i])**2)
    Dzs_fraction.append(Dznew)
    Dzs_rms_fraction.append(Dzs_stdv)
    Dzs_diff.append(Dzs_dyn[i]-Dparallel_dyn[i])
    Dzs_rms_diff.append(Dzs_stdv_diff)
    if spacing>1.9:
        spacings_plot.append(spacing)
        Dzs_fraction_plot.append(Dznew)
        Dzs_fraction_rms_plot.append(Dzs_stdv)
        ones.append(1)
    if spacing>2:
        spacings_few.append(spacing)
        Dzs_fraction_few.append(Dznew)
        Dzs_fraction_rms_few.append(Dzs_stdv)
    i+=1
    spacing = spacings_dyn[i]

ones = np.array(ones)
Dzs_diff = np.array(Dzs_diff)
Dzs_fraction = np.array(Dzs_fraction)
Dzs_rms_diff = np.array(Dzs_rms_diff)
spacings_few = np.array(spacings_few)
Dzs_fraction_few = np.array(Dzs_fraction_few)
Dzs_rms_fraction = np.array(Dzs_rms_fraction)
spacings_fraction = np.array(spacings_fraction)
Dzs_fraction_few_dyn = np.array(Dzs_fraction_plot)
Dzs_fraction_few_dyn_rms = np.array(Dzs_fraction_rms_plot)

Dzs_fraction_few_dyn_prms = Dzs_fraction_few_dyn+Dzs_fraction_few_dyn_rms
Dzs_fraction_few_dyn_mrms = Dzs_fraction_few_dyn-Dzs_fraction_few_dyn_rms


# Power law:
coeffs, covs = curve_fit(mypowerlaw,spacings_few, Dzs_fraction_few) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
Apl_dyn = coeffs[0]
bpl_dyn = coeffs[1]
rms_Apl_dyn = np.sqrt(covs[0,0])
rms_bpl_dyn = np.sqrt(covs[1,1])

fitpl_dyn   = mypowerlaw(spacings_fit,Apl_dyn,bpl_dyn)
fitpl_dyn_short = mypowerlaw(spacings_few,Apl_dyn,bpl_dyn)


#'''
fig, ax1 = plt.subplots(dpi=300)
ax1.fill_between(spacings_plot, Dzs_fraction_few_dyn_prms, Dzs_fraction_few_dyn_mrms,color='lightcoral',alpha=0.2)
ax1.plot(spacings_plot, Dzs_fraction_few_dyn, '-o',color='lightcoral',label='$D_{\perp, \mathregular{dyn}}/D_{\parallel, \mathregular{dyn}}$')
ax1.axhline(y=1,color='darkgray',linestyle=(0, (1, 1)))#,linewidth=1.0)
ax1.plot(spacings_fit,fitpl_dyn,'--',color='tab:blue',label=r'Power law fit')
ax1.set_xlabel(r'$d$ (nm)', fontsize=15)
ax1.set_ylabel(r'$D_{\perp}/D_{\parallel}$', fontsize=15)
ax1.legend(loc='upper center',ncol=3,fontsize=12)
plt.tight_layout()

ax2 = fig.add_axes([0.4,0.4,0.5,0.4])
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax2.loglog(spacings_plot, Dzs_fraction_few_dyn-1, 'o', color='lightcoral')
ax2.loglog(spacings_few, fitpl_dyn_short-1,'--', color='tab:blue')
ax2.set_xlabel(xlabel=r'$d$ (nm)',labelpad=-1.5,fontsize=12)
ax2.set_ylabel(ylabel=r'$D_\perp/D_\parallel$-1',labelpad=-3,fontsize=12)

plt.savefig(plotname_Dz_fraction_dyn)


### Stat, perp vs parallel

spacings_fraction = []
Dzs_fraction      = []
Dzs_diff          = []
Dzs_rms_fraction  = []
Dzs_rms_diff      = []
ones              = []
spacings_few      = []
Dzs_fraction_few  = []
Dzs_fraction_rms_few = []

i = 0
spacing = spacings_f[i]
while spacing<=25:
    spacings_fraction.append(spacing)    
    Dznew = Dzs_stat[i]/Dparallel_stat[i]
    # Ds, stdv
    Dzs_stdv = abs(Dznew)*np.sqrt((Dzs_stdv_stat[i]/Dzs_stat[i])**2+(Dparallel_stdv_stat[i]/Dparallel_stat[i])**2)
    Dzs_stdv_diff = np.sqrt((Dzs_stdv_stat[i])**2+(Dparallel_stdv_stat[i])**2)
    Dzs_fraction.append(Dznew)
    Dzs_rms_fraction.append(Dzs_stdv)
    Dzs_diff.append(Dzs_stat[i]-Dparallel_stat[i])
    Dzs_rms_diff.append(Dzs_stdv_diff)
    if spacing>1.5:
        spacings_few.append(spacing)
        Dzs_fraction_few.append(Dznew)
        Dzs_fraction_rms_few.append(Dzs_stdv)
        ones.append(1)
    i+=1
    spacing = spacings_stat[i]

ones = np.array(ones)
Dzs_diff = np.array(Dzs_diff)
spacings_few = np.array(spacings_few)
Dzs_fraction = np.array(Dzs_fraction)
Dzs_rms_diff = np.array(Dzs_rms_diff)
Dzs_fraction_few = np.array(Dzs_fraction_few)
Dzs_rms_fraction = np.array(Dzs_rms_fraction)
spacings_fraction = np.array(spacings_fraction)
Dzs_fraction_few_stat = np.array(Dzs_fraction_few)
Dzs_fraction_few_stat_rms = np.array(Dzs_fraction_rms_few)

Dzs_fraction_few_stat_prms = Dzs_fraction_few_stat+Dzs_fraction_few_stat_rms
Dzs_fraction_few_stat_mrms = Dzs_fraction_few_stat-Dzs_fraction_few_stat_rms

# Power law:
coeffs, covs = curve_fit(mypowerlaw,spacings_few, Dzs_fraction_few) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
Apl_stat = coeffs[0]
bpl_stat = coeffs[1]
rms_Apl_stat = np.sqrt(covs[0,0])
rms_bpl_stat = np.sqrt(covs[1,1])

fitpl_stat   = mypowerlaw(spacings_fit,Apl_stat,bpl_stat)
fitpl_stat_short = mypowerlaw(spacings_few,Apl_stat,bpl_stat)

fig, ax1 = plt.subplots(dpi=300)
ax1.fill_between(spacings_plot, Dzs_fraction_few_stat_prms, Dzs_fraction_few_stat_mrms,color='lightcoral',alpha=0.2)
ax1.plot(spacings_plot, Dzs_fraction_few_stat, '-o',color='lightcoral',label='$D_{\perp, \mathregular{stat}}/D_{\parallel, \mathregular{stat}}$')
ax1.axhline(y=1,color='darkgray',linestyle=(0, (1, 1)))#,linewidth=1.0)
ax1.plot(spacings_fit,fitpl_stat,'--',color='tab:blue',label=r'Power law fit')
ax1.set_xlabel(r'$d$ (nm)', fontsize=15)
ax1.set_ylabel(r'$D_{\perp}/D_{\parallel}$', fontsize=15)
ax1.legend(loc='lower center',ncol=3,fontsize=12)
ax1.axis([1.7,25.3,-1,11])
plt.tight_layout()

ax2 = fig.add_axes([0.4,0.5,0.4,0.4]) #fig.add_axes([0.28,0.5,0.4,0.4])
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax2.loglog(spacings_plot, Dzs_fraction_few_stat-1, 'o', color='lightcoral')
ax2.loglog(spacings_few,fitpl_stat_short-1,'--', color='tab:blue')
ax2.set_xlabel(xlabel=r'$d$ (nm)',labelpad=-1.5,fontsize=12)
ax2.set_ylabel(ylabel=r'$D_\perp/D_\parallel$-1',labelpad=-3,fontsize=12)


plt.savefig(plotname_Dz_fraction_stat)

### No stiff, perp vs parallel

spacings_fraction = []
Dzs_fraction      = []
Dzs_diff          = []
Dzs_rms_fraction  = []
Dzs_rms_diff      = []
ones              = []
spacings_few      = []
Dzs_fraction_few  = []
Dzs_fraction_rms_few = []

i = 0
spacing = spacings_nostiff[i]
while spacing<25:
    spacing = spacings_nostiff[i]
    spacings_fraction.append(spacing)   
    #print('i:',i, 'len(Dzs_nostiff):',len(Dzs_nostiff), 'len(Dparallel_nostiff):',len(Dparallel_nostiff)) 
    #print('Dzs_nostiff[',i,']:',Dzs_nostiff[i]) 
    #print('Dparallel_nostiff[',i,']:',Dparallel_nostiff[i]) 
    Dznew = Dzs_nostiff[i]/Dparallel_nostiff[i]
    # Ds, stdv
    Dzs_stdv = abs(Dznew)*np.sqrt((Dzs_stdv_nostiff[i]/Dzs_nostiff[i])**2+(Dparallel_stdv_nostiff[i]/Dparallel_nostiff[i])**2)
    Dzs_stdv_diff = np.sqrt((Dzs_stdv_nostiff[i])**2+(Dparallel_stdv_nostiff[i])**2)
    Dzs_fraction.append(Dznew)
    Dzs_rms_fraction.append(Dzs_stdv)
    Dzs_diff.append(Dzs_nostiff[i]-Dparallel_nostiff[i])
    Dzs_rms_diff.append(Dzs_stdv_diff)
    if spacing>1.5:
        spacings_few.append(spacing)
        Dzs_fraction_few.append(Dznew)
        Dzs_fraction_rms_few.append(Dzs_stdv)
        ones.append(1)
    i+=1

ones = np.array(ones)
Dzs_diff = np.array(Dzs_diff)
spacings_few = np.array(spacings_few)
Dzs_fraction = np.array(Dzs_fraction)
Dzs_rms_diff = np.array(Dzs_rms_diff)
Dzs_fraction_few = np.array(Dzs_fraction_few)
Dzs_rms_fraction = np.array(Dzs_rms_fraction)
spacings_fraction = np.array(spacings_fraction)
Dzs_fraction_few_nostiff = np.array(Dzs_fraction_few)
Dzs_fraction_few_nostiff_rms = np.array(Dzs_fraction_rms_few)

Dzs_fraction_few_nostiff_prms = Dzs_fraction_few_nostiff+Dzs_fraction_few_nostiff_rms
Dzs_fraction_few_nostiff_mrms = Dzs_fraction_few_nostiff-Dzs_fraction_few_nostiff_rms

# Power law:
coeffs, covs = curve_fit(mypowerlaw,spacings_few, Dzs_fraction_few)
Apl_nostiff = coeffs[0]
bpl_nostiff = coeffs[1]
rms_Apl_nostiff = np.sqrt(covs[0,0])
rms_bpl_nostiff = np.sqrt(covs[1,1])

fitpl_nostiff   = mypowerlaw(spacings_fit,Apl_nostiff,bpl_nostiff)
fitpl_nostiff_short = mypowerlaw(spacings_few,Apl_nostiff,bpl_nostiff)

fig, ax1 = plt.subplots(dpi=300)
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)
ax1.fill_between(spacings_few, Dzs_fraction_few_nostiff_prms, Dzs_fraction_few_nostiff_mrms,color='lightcoral',alpha=0.2)
ax1.plot(spacings_few, Dzs_fraction_few_nostiff, '-o',color='lightcoral',label='$D_{\perp, \mathdefault{no stiff.}}/D_{\parallel, \mathdefault{no stiff.}}$')
ax1.axhline(y=1,color='darkgray',linestyle=(0, (1, 1)))
ax1.plot(spacings_fit,fitpl_nostiff,'--',color='tab:blue',label=r'Power law fit')
ax1.set_xlabel(r'$d$ (nm)', fontsize=15)
ax1.set_ylabel(r'$D_{\perp}/D_{\parallel}$', fontsize=15)
ax1.legend(loc='upper center',ncol=2,fontsize=12)
plt.tight_layout()

ax2 = fig.add_axes([0.4,0.4,0.5,0.4])
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax2.loglog(spacings_few, Dzs_fraction_few_nostiff-1, 'o', color='lightcoral',label='$D_{\perp, \mathregular{no stiff.}}/D_{\parallel, \mathregular{no stiff.}}$')
ax2.loglog(spacings_few,fitpl_nostiff_short-1,'--', color='tab:blue',label=r'Power law fit')
ax2.set_xlabel(xlabel=r'$d$ (nm)',labelpad=-1.5,fontsize=12)
ax2.set_ylabel(ylabel=r'$D_\perp/D_\parallel$-1',labelpad=-3,fontsize=12)


plt.savefig(plotname_Dz_fraction_nostiff)

### Straight, perp vs parallel
spacings_fraction = []
Dzs_fraction      = []
Dzs_diff          = []
Dzs_rms_fraction  = []
Dzs_rms_diff      = []
ones              = []
spacings_few      = []
Dzs_fraction_few  = []
Dzs_fraction_rms_few = []

spacings_plot     = []
Dzs_fraction_plot = []
Dzs_fraction_rms_plot = []

i = 0
spacing = spacings_f[i]
while spacing<25:
    spacing = spacings_f[i]
    spacings_fraction.append(spacing)   
    Dznew = Dzs_f[i]/Dparallel_f[i]
    # Ds, stdv
    Dzs_stdv = abs(Dznew)*np.sqrt((Dzs_stdv_f[i]/Dzs_f[i])**2+(Dparallel_stdv_f[i]/Dparallel_f[i])**2)
    Dzs_stdv_diff = np.sqrt((Dzs_stdv_f[i])**2+(Dparallel_stdv_f[i])**2)
    Dzs_fraction.append(Dznew)
    Dzs_rms_fraction.append(Dzs_stdv)
    Dzs_diff.append(Dzs_f[i]-Dparallel_f[i])
    Dzs_rms_diff.append(Dzs_stdv_diff)
    if spacing>1.9:
        spacings_plot.append(spacing)
        Dzs_fraction_plot.append(Dznew)
        Dzs_fraction_rms_plot.append(Dzs_stdv)
        ones.append(1)
    if spacing>3:
        spacings_few.append(spacing)
        Dzs_fraction_few.append(Dznew)
        Dzs_fraction_rms_few.append(Dzs_stdv)
    i+=1

Dzs_diff = np.array(Dzs_diff)
spacings_few = np.array(spacings_few)
Dzs_fraction = np.array(Dzs_fraction)
Dzs_rms_diff = np.array(Dzs_rms_diff)
Dzs_fraction_few = np.array(Dzs_fraction_few)
Dzs_rms_fraction = np.array(Dzs_rms_fraction)
spacings_fraction = np.array(spacings_fraction)
Dzs_fraction_few_straight = np.array(Dzs_fraction_plot)
Dzs_fraction_few_straight_rms = np.array(Dzs_fraction_rms_plot)

Dzs_fraction_few_straight_prms = Dzs_fraction_few_straight+Dzs_fraction_few_straight_rms
Dzs_fraction_few_straight_mrms = Dzs_fraction_few_straight-Dzs_fraction_few_straight_rms

# Power law:
coeffs, covs = curve_fit(mypowerlaw,spacings_few, Dzs_fraction_few) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
Apl_f = coeffs[0]
bpl_f = coeffs[1]
rms_Apl_f = np.sqrt(covs[0,0])
rms_bpl_f = np.sqrt(covs[1,1])

fitpl_f   = mypowerlaw(spacings_fit,Apl_f,bpl_f)
fitpl_f_short = mypowerlaw(spacings_few,Apl_f,bpl_f)

fig, ax1 = plt.subplots(dpi=300)
ax1.fill_between(spacings_plot, Dzs_fraction_few_straight_prms, Dzs_fraction_few_straight_mrms,color='lightcoral',alpha=0.2)
ax1.plot(spacings_plot, Dzs_fraction_few_straight, '-o',color='lightcoral',label='$D_{\perp, \mathregular{straight}}/D_{\parallel, \mathregular{straight}}$')
ax1.axhline(y=1,color='darkgray',linestyle=(0, (1, 1)))
ax1.plot(spacings_fit,fitpl_f,'--',color='tab:blue',label=r'Power law fit')
ax1.set_xlabel(r'$d$ (nm)', fontsize=15)
ax1.set_ylabel(r'$D_{\perp}/D_{\parallel}$', fontsize=15)
ax1.legend(loc='lower center',ncol=2,fontsize=12)
plt.tight_layout()

ax2 = fig.add_axes([0.5,0.5,0.4,0.4])#([0.35, 0.35, 0.45, 0.45])
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax2.loglog(spacings_plot, Dzs_fraction_few_straight-1, 'o', color='lightcoral')
ax2.loglog(spacings_few,fitpl_f_short-1,'--', color='tab:blue')
ax2.set_xlabel(xlabel=r'$d$ (nm)',labelpad=-1.5,fontsize=12)#14)
ax2.set_ylabel(ylabel=r'$D_\perp/D_\parallel$-1',labelpad=-3,fontsize=12)#14)

plt.savefig(plotname_Dz_fraction_forest)


print('########################## POWER LAW ###########################################\n')
print("Apl_dyn:",Apl_dyn, ", bpl_dyn:",bpl_dyn)
print("rms_Apl_dyn:",rms_Apl_dyn, ", rms_bpl_dyn:",rms_bpl_dyn)
print("----------------------------------")
print("Apl_stat:",Apl_stat, ", bpl_stat:",bpl_stat)
print("rms_Apl_stat:",rms_Apl_stat, ", rms_bpl_stat:",rms_bpl_stat)
print("----------------------------------")
print("Apl_nostiff:",Apl_nostiff, ", bpl_nostiff:",bpl_nostiff)
print("rms_Apl_nostiff:",rms_Apl_nostiff, ", rms_bpl_nostiff:",rms_bpl_nostiff)
print("----------------------------------")
print("Apl_f:",Apl_f, ", bpl_f:",bpl_f)
print("rms_Apl_f:",rms_Apl_f, ", rms_bpl_f:",rms_bpl_f)



plt.show()
