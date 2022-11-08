import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob
import copy

Nrms = 10 # Number of averages 

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

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    print('avgx:',avgx)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

#
damp = 10
psigma       = 3
confignrs    = np.arange(1,1001)
Nsteps       = 2001 
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime

filenames = []
ds = np.array([4,6,8,10]) # These can possibly vary from sigma to sigma

for d in ds:
    pathstring = 'Diffusion_bead_near_grid/Spacing'+str(d)+'/damp10_diffseedLgv/Brush/Sigma_bead_'+str(psigma)+'/Nocut/diffusionconfig1to1000_better_rms_Nestimates10.txt'
    filenames.append(pathstring)


N = len(ds)

DR   = np.zeros(N)
Dz   = np.zeros(N)
Dpar = np.zeros(N)
DR_rms   = np.zeros(N)
Dz_rms   = np.zeros(N)
Dpar_rms = np.zeros(N)

print('filenames')

outfilename = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_3/Nocut/D_vs_d_sigma'+str(psigma)+'.txt'
outfile = open(outfilename, 'w')

i = 0
for filename in filenames:
    infile = open(filename,'r')
    header = infile.readline()
    data   = infile.readline()
    words  = data.split()
    
    #    [0]     [1]    [2]     [3]      [4]      [5]
    # (DR_avg, DR_rms, Dz_avg, Dz_rms, Dpar_avg, Dpar_rms)
    thisDz = float(words[2])
    thisDz_rms = float(words[3])
    outfile.write('%i %.5e %.5e\n' % (ds[i],thisDz,thisDz_rms))
    print('ds[',i,']:',ds[i])
    print('thisDz:',thisDz)
    print('thisDz_rms:',thisDz_rms)
    
    DR[i] = float(words[0])
    Dz[i] = thisDz
    Dpar[i] = float(words[4])
    DR_rms[i] = float(words[1])
    Dz_rms[i] = thisDz_rms
    Dpar_rms[i] = float(words[5])
    
    infile.close()
    i+=1

outfile.close()
