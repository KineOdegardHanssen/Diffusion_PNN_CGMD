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

# change the default font family
plt.rcParams.update({'font.family':'Arial'})
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

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
#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
long    = False
spacing = 1.5
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)

# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
zhigh          = 250
zlow           = -50
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)
maxz_av      = 0
filescounter = 0
if long==True:
    Nsteps   = 10001
else:
    Nsteps   = 2001
unitlength   = 1e-9     # m
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime
Ninch = 101
Nch   = 9 

endlocation_in           = 'Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/'
endlocation              = endlocation_in +'Nocut/'
filestext                = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
# Text files
outfilename              = endlocation+'corrlen_'+filestext+'_nocut.txt'

skippedfiles = 0

xsiaverage  = 0
xsirmsarray = []
for confignr in confignrs:
    print('On config number:', confignr)
    if long==True:
        infilename_all  = endlocation_in+'long/'+'all_confignr'+str(confignr)+'_long.lammpstrj'
    else:
        infilename_all  = endlocation_in+'all_confignr'+str(confignr)+'.lammpstrj'
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    try:
        infile_all = open(infilename_all, "r")
    except:
        try:
            infilename_all = endlocation_in+'long/'+'all_confignr'+str(confignr)+'.lammpstrj'
            infile_all = open(infilename_all, "r")
        except:
            print('Oh, lammpstrj-file! Where art thou?')
            skippedfiles += 1
            continue # Skipping this file if it does not exist
    # Moving on, if the file
    filescounter += 1
    lines = infile_all.readlines()
    # Getting the number of lines, etc.
    totlines = len(lines)         # Total number of lines
    lineend = totlines-1          # Index of last element
    
    # Extracting the number of atoms:
    words = lines[3].split()
    Nall = int(words[0])
    N    = Nall
    
    skiplines   = 9             # If we hit 'ITEM:', skip this many steps
    skipelem    = 0
    sampleevery = 0
    i           = int(math.ceil(skipelem*(Nall+9)))
    skiplines  += (Nall+skiplines)*sampleevery
    
    maxz = -1000 # No z-position is this small.
    
    positions = np.zeros((Nch,Ninch,3))
    
    time_start = time.process_time()
    counter = 0
    while i<totlines:
        words = lines[i].split()
        if (words[0]=='ITEM:'):
            if words[1]=='TIMESTEP':
                words2 = lines[i+1].split() # The time step is on the next line
                t = float(words2[0])
                counters = []
                for j in range(Nch):
                    counters.append(0)
                i+=skiplines
            elif words[1]=='NUMBER': # Safeguards, should not kick in
                i+=7
            elif words[1]=='BOX':
                i+=5
            elif words[1]=='ATOMS':
                i+=1
        elif len(words)<8:
            i+=1
        else:
            if t>0:
                # Find properties
                # Order:  id  type mol ux  uy  uz  vx  vy   vz
                #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
                ind      = int(words[0])-1 # Atom ids go from zero to N-1.
                atomtype = int(words[1]) 
                if atomtype==1 or atomtype==2: # Chain bead
                    x = float(words[3])
                    y = float(words[4])
                    z = float(words[5])
                    molID = int(words[2])-1 # Different indexation
                    positions[molID,counters[molID],0]=x
                    positions[molID,counters[molID],1]=y
                    positions[molID,counters[molID],2]=z
                    counters[molID] += 1
            i+=1
    infile_all.close()
    
    # Loop over all chains and beads to find the closest one for each bead
    for i in range(Nch): # For each chain
        for j in range(Ninch): # For each bead
            smallestdist = 1e8        # Find the closest bead on another chain  (No distance is this big)
            thispos = positions[i,j,:]
            for k in range(Nch): 
                if i!=k:                # Loop over all the other chains
                    for l in range(Ninch): # Loop over all beads in that chain
                        candpos    = positions[k,l,:]
                        difference = thispos-candpos
                        distance2  = np.dot(difference,difference)
                        if distance2<smallestdist:  # Get the smallest distance
                            smallestdist = distance2
            # Accumulate average:
            smallestdist = np.sqrt(smallestdist)
            xsirmsarray.append(smallestdist)     # Store the shortest distance to a bead on another chain. Do this for all beads
   
xsirmsarray = np.array(xsirmsarray)
print('max(xsirmsarray):',max(xsirmsarray))
print('min(xsirmsarray):',min(xsirmsarray))
print('np.mean(xsirmsarray):',np.mean(xsirmsarray))
print('xsirmsarray:',xsirmsarray)
Nxsis = len(xsirmsarray)
xsiaverage = np.mean(xsirmsarray)

xsirms = 0
for i in range(Nxsis):
    xsirms = (xsiaverage-xsirmsarray[i])**2
xsirms = np.sqrt(xsirms/(Nxsis-1))

outfile = open(outfilename,'w')
outfile.write('%.5e %.5e' % (xsiaverage,xsirms))
outfile.close()
