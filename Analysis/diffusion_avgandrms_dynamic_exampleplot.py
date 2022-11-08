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

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

plottest = False # True

Nrms = 10 # Number of averages 

DR_estimates = np.zeros(Nrms)
Dz_estimates = np.zeros(Nrms)
Dpar_estimates = np.zeros(Nrms)
# Constant terms (for plotting purposes)
bR_estimates = np.zeros(Nrms)
bz_estimates = np.zeros(Nrms)
bpar_estimates = np.zeros(Nrms)

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
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

#
damp = 10
# Input parameters for file selection:
popup_plots = False
long    = False
spacing = 5
psigma  = 1

# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
zhigh          = 250
zlow           = -50
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)
N_per_rms    = int(Nseeds/Nrms)
maxz_av      = 0
filescounter = 0
if long==True:
    Nsteps   = 10001
else:
    Nsteps   = 2001
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime

endlocation_in = 'Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/'
endlocation    = endlocation_in +'Nocut/'
filestext      = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

# Plots
plotname       = endlocation+filestext+'_nocut_better_rms_Nestimates%i_procedure.png' % Nrms


## Getting interval to perform the fit:
rangefilename = endlocation+'indices_for_fit_better_rms.txt'
rangefile     = open(rangefilename,'r')
lines         = rangefile.readlines()  
startindex_R  = int(lines[0].split()[1])
endindex_R    = int(lines[0].split()[2])
startindex_z  = int(lines[1].split()[1])
endindex_z    = int(lines[1].split()[2])
startindex_p  = int(lines[2].split()[1])
endindex_p    = int(lines[2].split()[2])
rangefile.close()


### Setting arrays
# These are all squared:
# All together:
alldxs     = []
alldys     = []
alldzs     = []
alldxs_cut = []
alldys_cut = []
alldzs_cut = []
## Averages:
# Nocut
tot_averageR2s            = np.zeros(Nsteps) # Distances, squared
tot_averagedx2s           = np.zeros(Nsteps)
tot_averagedy2s           = np.zeros(Nsteps)
tot_averagedz2s           = np.zeros(Nsteps)
tot_averagedparallel2     = np.zeros(Nsteps) # Distances, squared
tot_average_counter       = np.zeros(Nsteps)
# Cut
tot_averageR2s_cut        = np.zeros(Nsteps) # Distances, squared
tot_averagedx2s_cut       = np.zeros(Nsteps)
tot_averagedy2s_cut       = np.zeros(Nsteps)
tot_averagedz2s_cut       = np.zeros(Nsteps)
tot_averagedparallel2_cut = np.zeros(Nsteps) # Distances, squared
tot_average_counter_cut   = np.zeros(Nsteps)

# Separated by seed:
Rs_byseed    = []
dxs_byseed   = []
dys_byseed   = []
dzs_byseed   = []
gamma_byseed = []
times_byseed = []
Nins         = []
#
# This is not squared, obviously:
alltimes  = []

for j in range(Nrms):
    configbase = N_per_rms*j+1
    # Nocut
    averageR2s            = np.zeros(Nsteps) # Distances, squared
    averagedx2s           = np.zeros(Nsteps)
    averagedy2s           = np.zeros(Nsteps)
    averagedz2s           = np.zeros(Nsteps)
    averagedxs            = np.zeros(Nsteps) # Distances, non-squared
    averagedys            = np.zeros(Nsteps)
    averagedzs            = np.zeros(Nsteps)
    averagedparallel2     = np.zeros(Nsteps) # Distances, squared
    average_counter       = np.zeros(Nsteps)
    # Cut
    averageR2s_cut        = np.zeros(Nsteps) # Distances, squared
    averagedx2s_cut       = np.zeros(Nsteps)
    averagedy2s_cut       = np.zeros(Nsteps)
    averagedz2s_cut       = np.zeros(Nsteps)
    averagedxs_cut        = np.zeros(Nsteps) # Distances, non-squared
    averagedys_cut        = np.zeros(Nsteps)
    averagedzs_cut        = np.zeros(Nsteps)
    averagedparallel2_cut = np.zeros(Nsteps) # Distances, squared
    average_counter_cut   = np.zeros(Nsteps)
    for k in range(N_per_rms):
        confignr = k+configbase
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
                infilename_all = endlocation+'long/'+'all_confignr'+str(confignr)+'.lammpstrj' 
                infile_all = open(infilename_all, "r")
            except:
                print('Oh, lammpstrj-file! Where art thou?')
                continue # Skipping this file if it does not exist
        # Moving on, if the file
        lines = infile_all.readlines()
        # Getting the number of lines, etc.
        totlines = len(lines)         # Total number of lines
        lineend = totlines-1          # Index of last element
        
        # Extracting the number of atoms:
        words = lines[3].split()
        Nall = int(words[0])
        N    = Nall
        
        skiplines   = 9             # If we hit 'ITEM:', skip this many steps.
        skipelem    = 0
        sampleevery = 0
        i           = int(math.ceil(skipelem*(Nall+9)))
        skiplines  += (Nall+skiplines)*sampleevery
        
        maxz = -1000 # No z-position is this small.
        
        time_start = time.process_time()
        counter = 0
        while i<totlines:
            words = lines[i].split()
            if words[0]=='ITEM:':
                if words[1]=='TIMESTEP':
                    i+=skiplines
                elif words[1]=='NUMBER': # Safeguard, these should never kick in. 
                    i+=7
                elif words[1]=='BOX':
                    i+=5
                elif words[1]=='ATOMS':
                    i+=1
            elif len(words)<8:
                i+=1
            else:
                # Find properties
                # Order:  id  type mol ux  uy  uz  vx  vy   vz
                #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
                ind      = int(words[0])-1 # Atom ids go from zero to N-1.
                atomtype = int(words[1]) 
                #molID    = int(words[2])
                z        = float(words[5])
                if atomtype==2: # Moving polymer bead. Test if this is larger than maxz:
                    if z>maxz:
                        maxz = z
                i+=1
        extent_polymers = maxz
        maxz_av += maxz
        infile_all.close()
        
        if long==True:
            infilename_free = endlocation_in+'long/'+'freeatom_confignr'+str(confignr)+'_long.lammpstrj'
        else:
            infilename_free = endlocation_in+'freeatom_confignr'+str(confignr)+'.lammpstrj'
    
        # Read in:
        #### Automatic part
        ## Find the extent of the polymers: Max z-coord of beads in the chains
        try:
            infile_free = open(infilename_free, "r")
        except:
            try:
                infilename_free = endlocation_in+'long/'+'freeatom_confignr'+str(confignr)+'.lammpstrj'
                infile_free = open(infilename_free, "r")
            except:
                print('Oh, lammpstrj-file! Where art thou?')
                continue # Skipping this file if it does not exist
        # Moving on, if the file
        
        ## Find the position of the free bead:
        lines = infile_free.readlines()
        # Getting the number of lines, etc.
        totlines = len(lines)         # Total number of lines
        lineend = totlines-1          # Index of last element
        
        # Extracting the number of atoms:
        words = lines[3].split()
        Nall = int(words[0])
        N    = Nall
        
        skiplines   = 9             # If we hit 'ITEM:', skip this many steps.
        skipelem    = 0
        sampleevery = 0
        i           = int(math.ceil(skipelem*(Nall+9)))
        skiplines  += (Nall+skiplines)*sampleevery
        
        # Setting arrays for treatment:
        
        # Setting arrays for treatment:
        positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation.
        times     = [] 
        
        maxz = -1000 # No z-position is this small.
        
        time_start = time.process_time()
        counter = 0
        zs_fortesting = []
        while i<totlines:
            words = lines[i].split()
            if (words[0]=='ITEM:'):
                if words[1]=='TIMESTEP':
                    words2 = lines[i+1].split() # The time step is on the next line
                    t = float(words2[0])
                    times.append(t)
                    i+=skiplines
                elif words[1]=='NUMBER': # These will never kick in. 
                    i+=1
                elif words[1]=='BOX':
                    i+=1
                elif words[1]=='ATOMS':
                    i+=1
            elif len(words)<8:
                i+=1
            else:
                # Find properties
                # Order:  id  type mol ux  uy  uz  vx  vy   vz
                #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
                x        = float(words[3])
                y        = float(words[4])
                z        = float(words[5])
                positions.append(np.array([x,y,z]))
                zs_fortesting.append(z)
                counter+=1
                i+=1
        infile_free.close()
        if len(zs_fortesting)<Nsteps:
            continue
        dt = (times[1]-times[0])*timestepsize
        times_single = np.arange(Nsteps)
        times_single_real = np.arange(Nsteps)*dt
        
        time_end = time.process_time()
        
        if max(zs_fortesting)>zhigh:
            continue
        if min(zs_fortesting)<zlow:
            continue
        
        filescounter += 1
        
        Nin   = 0
        maxzpol = 0
        pos_inpolymer = []
        for i in range(counter):
            thesepos = positions[i]
            z        = thesepos[2]
            if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
                break
            else:
                pos_inpolymer.append(thesepos)
                if z>maxzpol:
                    maxzpol = z
                Nin+=1
        #######
        R_temp  = []
        dx_temp = []
        dy_temp = []
        dz_temp = []
        step_temp = []
        
        R_temp.append(0)
        dx_temp.append(0)
        dy_temp.append(0)
        dz_temp.append(0)
        step_temp.append(0)
        startpos = positions[0]
        for i in range(1,Nsteps):       
            this = positions[i]
            dist = this-startpos
            dx   = dist[0]         # Signed
            dy   = dist[1]
            dz   = dist[2]
            R2   = np.dot(dist,dist) # Squared
            dx2  = dx*dx
            dy2  = dy*dy
            dz2  = dz*dz
            # Averages, this batch of files:
            averageR2s[i] +=R2   # Distance
            averagedx2s[i]+=dx2
            averagedy2s[i]+=dy2
            averagedz2s[i]+=dz2
            averagedxs[i] +=dx   # Distance, signed
            averagedys[i] +=dy
            averagedzs[i] +=dz
            averagedparallel2[i] += dx2+dy2 # Distance
            average_counter[i] +=1
            # Total averages: ############### Only need these for plotting
            tot_averageR2s[i] +=R2   # Distances, squared
            tot_averagedx2s[i] +=dx2
            tot_averagedy2s[i] +=dy2
            tot_averagedz2s[i] +=dz2
            tot_averagedparallel2[i] += dx2+dy2 # Distances, squared
            tot_average_counter[i]+=1
            ################################
            # Separated by seed:
            R_temp.append(R2)
            dx_temp.append(dx2)
            dy_temp.append(dy2)
            dz_temp.append(dz2)
            step_temp.append(i)
        
        for i in range(1,Nin):       
            this = pos_inpolymer[i]
            dist = this-startpos
            dx   = dist[0]         # Signed
            dy   = dist[1]
            dz   = dist[2]
            R2   = np.dot(dist,dist) # Squared
            dx2  = dx*dx
            dy2  = dy*dy
            dz2  = dz*dz
            # Averages, this batch of files:
            averageR2s_cut[i] +=R2   # Distance
            averagedx2s_cut[i]+=dx2
            averagedy2s_cut[i]+=dy2
            averagedz2s_cut[i]+=dz2
            averagedxs_cut[i] +=dx   # Distance, signed
            averagedys_cut[i] +=dy
            averagedzs_cut[i] +=dz
            averagedparallel2_cut[i] += dx2+dy2 # Distance
            average_counter_cut[i] +=1
            # Total averages: ############### For plotting
            tot_averageR2s_cut[i] +=R2   # Distances, squared
            tot_averagedx2s_cut[i] +=dx2
            tot_averagedy2s_cut[i] +=dy2
            tot_averagedz2s_cut[i] +=dz2
            tot_averagedparallel2_cut[i] += dx2+dy2 # Distances, squared
            tot_average_counter_cut[i]+=1
        
        Nins.append(Nin)
    
    # Calculating averages:
    
    Ninbrush = Nsteps
    for i in range(1,Nsteps):
        counter = average_counter[i]
        if counter!=0:
            # Distance squared
            averageR2s[i]/=counter
            averagedx2s[i]/=counter
            averagedy2s[i]/=counter
            averagedz2s[i]/=counter
            averagedparallel2[i]/= counter
            # Distance, signed:
            averagedxs[i]/=counter
            averagedys[i]/=counter
            averagedzs[i]/=counter
        else:
           Ninbrush = i-1
           break
    Ninbrush_cut = Nin
    for i in range(1,Nin):
        counter = average_counter_cut[i]
        if counter!=0:
            # Distance squared
            averageR2s_cut[i]/=counter
            averagedx2s_cut[i]/=counter
            averagedy2s_cut[i]/=counter
            averagedz2s_cut[i]/=counter
            averagedparallel2_cut[i]/= counter
            # Distance, signed:
            averagedxs_cut[i]/=counter
            averagedys_cut[i]/=counter
            averagedzs_cut[i]/=counter
        else:
           Ninbrush_cut = i-1
           break
    
    times_single      = times_single[0:Ninbrush]
    times_single_real = times_single_real[0:Ninbrush]
    # Distance squared
    averageR2s  = averageR2s[0:Ninbrush]
    averagedx2s = averagedx2s[0:Ninbrush]
    averagedy2s = averagedy2s[0:Ninbrush]
    averagedz2s = averagedz2s[0:Ninbrush]
    averagedparallel2 = averagedparallel2[0:Ninbrush]
    averageR2s_cut  = averageR2s_cut[0:Ninbrush_cut]
    averagedx2s_cut = averagedx2s_cut[0:Ninbrush_cut]
    averagedy2s_cut = averagedy2s_cut[0:Ninbrush_cut]
    averagedz2s_cut = averagedz2s_cut[0:Ninbrush_cut]
    averagedparallel2_cut = averagedparallel2_cut[0:Ninbrush_cut]
    # Distance to the first power
    averagedxs = averagedxs[0:Ninbrush]
    averagedys = averagedys[0:Ninbrush]
    averagedzs = averagedzs[0:Ninbrush]
    
    ### Average, SI units:
    ## Distance:
    # Nocut
    averageR2s_SI  = averageR2s*unitlength**2
    averagedx2s_SI = averagedx2s*unitlength**2
    averagedy2s_SI = averagedy2s*unitlength**2
    averagedz2s_SI = averagedz2s*unitlength**2
    averagedparallel2_SI = averagedparallel2*unitlength**2
    # Cut
    averageR2s_cut_SI  = averageR2s_cut*unitlength**2
    averagedx2s_cut_SI = averagedx2s_cut*unitlength**2
    averagedy2s_cut_SI = averagedy2s_cut*unitlength**2
    averagedz2s_cut_SI = averagedz2s_cut*unitlength**2
    averagedparallel2_cut_SI = averagedparallel2_cut*unitlength**2
    
    # Insert interval
    coeffs_poly, covs = polyfit(times_single_real[startindex_R:endindex_R], averageR2s_SI[startindex_R:endindex_R], 1, full=False, cov=True) # Using all the data in the fit 
    a_poly_SI = coeffs_poly[0]
    b_poly_SI = coeffs_poly[1]
    rms_D_poly_SI = np.sqrt(covs[0,0])/6.
    rms_b_poly_SI = np.sqrt(covs[1,1])
    D_poly_SI = a_poly_SI/6.
    
    DR_estimates[j] = D_poly_SI
    bR_estimates[j] = b_poly_SI
    
    fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI
    if plottest==True:
       plt.figure(figsize=(6,5))
       plt.plot(times_single_real, averageR2s_SI)
       plt.plot(times_single_real, fit_poly_SI, '--')
       plt.xlabel('Time [s]')
       plt.ylabel('$<dR^2>$')
       plt.title('Test, R')
       plt.show()
    
    coeffs_poly, covs = polyfit(times_single_real[startindex_p:endindex_p], averagedparallel2_SI[startindex_p:endindex_p], 1, full=False, cov=True) # Using all the data in the fit 
    a_poly_SI = coeffs_poly[0]
    b_poly_SI = coeffs_poly[1]
    rms_D_poly_SI = np.sqrt(covs[0,0])/6.
    rms_b_poly_SI = np.sqrt(covs[1,1])
    D_poly_SI = a_poly_SI/4.
    
    Dpar_estimates[j] = D_poly_SI
    bpar_estimates[j] = b_poly_SI
    
    fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI
    if plottest==True:
       plt.figure(figsize=(6,5))
       plt.plot(times_single_real, averagedparallel2_SI)
       plt.plot(times_single_real, fit_poly_SI, '--')
       plt.xlabel('Time [s]')
       plt.ylabel('$<dparallel^2>$')
       plt.title('Test, parallel')
       plt.show()
    
    coeffs_poly, covs = polyfit(times_single_real[startindex_z:endindex_z], averagedz2s_SI[startindex_z:endindex_z], 1, full=False, cov=True) # Using all the data in the fit 
    a_poly_SI = coeffs_poly[0]
    b_poly_SI = coeffs_poly[1]
    rms_D_poly_SI = np.sqrt(covs[0,0])/6.
    rms_b_poly_SI = np.sqrt(covs[1,1])
    D_poly_SI = a_poly_SI/2.
    
    Dz_estimates[j] = D_poly_SI
    bz_estimates[j] = b_poly_SI
    
    fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI
    if plottest==True:
       plt.figure(figsize=(6,5))
       plt.plot(times_single_real, averagedz2s_SI)
       plt.plot(times_single_real, fit_poly_SI, '--')
       plt.xlabel('Time [s]')
       plt.ylabel('$<dz^2>$')
       plt.title('Test, z')
       plt.show()
    
for i in range(1,Nsteps):
    counter = tot_average_counter[i]
    if counter!=0:
        # Distance squared
        tot_averageR2s[i]/=counter
        tot_averagedx2s[i]/=counter
        tot_averagedy2s[i]/=counter
        tot_averagedz2s[i]/=counter
        tot_averagedparallel2[i]/= counter
    else:
        Ninbrush = i-1
        break
for i in range(1,Nin):
    counter = tot_average_counter_cut[i]
    if counter!=0:
        # Distance squared
        tot_averageR2s_cut[i]/=counter
        tot_averagedx2s_cut[i]/=counter
        tot_averagedy2s_cut[i]/=counter
        tot_averagedz2s_cut[i]/=counter
        tot_averagedparallel2_cut[i]/= counter
    else:
        Ninbrush_cut = i-1
        break

tot_averageR2s_SI  = tot_averageR2s*unitlength**2
tot_averagedx2s_SI = tot_averagedx2s*unitlength**2
tot_averagedy2s_SI = tot_averagedy2s*unitlength**2
tot_averagedz2s_SI = tot_averagedz2s*unitlength**2
tot_averagedparallel2_SI = tot_averagedparallel2*unitlength**2

tot_averageR2s_cut_SI  = tot_averageR2s_cut*unitlength**2
tot_averagedx2s_cut_SI = tot_averagedx2s_cut*unitlength**2
tot_averagedy2s_cut_SI = tot_averagedy2s_cut*unitlength**2
tot_averagedz2s_cut_SI = tot_averagedz2s_cut*unitlength**2
tot_averagedparallel2_cut_SI = tot_averagedparallel2_cut*unitlength**2


DR_avg, DR_rms = avg_and_rms(DR_estimates)
Dz_avg, Dz_rms = avg_and_rms(Dz_estimates)
Dpar_avg, Dpar_rms = avg_and_rms(Dpar_estimates)
# Constant term for plotting purposes:
bR_avg, bR_rms = avg_and_rms(bR_estimates)
bz_avg, bz_rms = avg_and_rms(bz_estimates)
bpar_avg, bpar_rms = avg_and_rms(bpar_estimates)

fit_DR_SI   = 6*DR_avg*times_single_real+bR_avg
fit_Dz_SI   = 2*Dz_avg*times_single_real+bz_avg
fit_Dpar_SI = 4*Dpar_avg*times_single_real+bpar_avg

# Plot
fig, ax = plt.subplots(figsize=(8,5),dpi=300)
# Cut
ax.plot(times_single_real, tot_averagedz2s_SI, label=r'Data, no cut')
ax.plot(times_single_real, tot_averagedz2s_cut_SI, label=r'Data, cut')
# No cut
ax.plot(times_single_real, fit_Dz_SI, '--', label=r'Fit')
ax.set_xlabel(r'Time (s)')
ax.set_ylabel(r'$<dz^2>$ (m$^2$)')
ax.set_title(r'$<dz^2>$ in brush, d = %.2f$\sigma$' % spacing)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
plt.legend(loc='lower right')

ax2 = fig.add_axes([0.2,0.55,0.3,0.27])
# Cut
# No cut
ax2.plot(times_single_real[startindex_R:endindex_R], tot_averagedz2s_SI[startindex_R:endindex_R], label=r'No cut')
ax2.plot(times_single_real[startindex_R:endindex_R], tot_averagedz2s_cut_SI[startindex_R:endindex_R], label=r'Cut')
ax2.plot(times_single_real[startindex_R:endindex_R], fit_Dz_SI[startindex_R:endindex_R], '--', label=r'Fit')
ax2.set_xlabel(r'Time (s)')
ax2.set_ylabel(r'$<dz^2>$ (m$^2$)')
plt.savefig(plotname)


plt.show()

print('spacing:', spacing)
print('psigma:', psigma)
