import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plot_in_tiff = True

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

def theory(a):
    k = 1.38e-23 # J/K
    T = 310      # K
    kT = k*T
    unittime = 2.38e-11  # s, unit time step
    a0 = 5e-10           # a at sigma = 1
    mass = 8.1e-25       # kg, unit mass
    damp = 10*unittime   # damp in SI units
    eta = mass/(6*np.pi*damp*a0)
    gamma = 6*np.pi*eta*a*1e-9
    return kT/gamma

Nintervals = 10

# Input parameters for file selection:
psigmas  = [0.25,0.5,1,1.5,2,3,5,10]
spacing  = 10
damp     = 10
N        = len(psigmas)
psigmas_base  = psigmas

# Input booleans for file selection:
bulkdiffusion = True
substrate     = False
confignrs     = np.arange(1,1001)

# Ds
DRs = np.zeros(N)
Dxs = np.zeros(N)
Dys = np.zeros(N)
Dzs = np.zeros(N)
Dparallel = np.zeros(N)
# Ds, stdv
DRs_stdv = np.zeros(N)
Dxs_stdv = np.zeros(N)
Dys_stdv = np.zeros(N)
Dzs_stdv = np.zeros(N)
Dparallel_stdv = np.zeros(N)

# bs
bRs = np.zeros(N)
bxs = np.zeros(N)
bys = np.zeros(N)
bzs = np.zeros(N)
bparallel = np.zeros(N)
# bs, stdv
bRs_stdv = np.zeros(N)
bxs_stdv = np.zeros(N)
bys_stdv = np.zeros(N)
bzs_stdv = np.zeros(N)
bparallel_stdv = np.zeros(N)


if bulkdiffusion==True:
    parentfolder = 'Bulk/'
    filestext    = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
    systemtype   = 'bulk'
    if substrate==True:
        parentfolder = 'Bulk_substrate/'
        systemtype   = 'substrate'
else:
    parentfolder = 'Brush/'
    systemtype   = 'brush'
    filestext    = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation_out   = 'Diffusion_bead_near_grid/D_vs_d/'+parentfolder+'Varysigmas/'
outfilename       = endlocation_out+'D_vs_sigma_better_rms_Nestimates%i_pmass1.txt' % Nintervals
plotname          = endlocation_out+'D_vs_sigma_bulk_better_rms_Nestimates%i_pmass1.png' % Nintervals
plotname_two      = endlocation_out+'D_vs_sigma_bulk_better_rms_Nestimates%i_twoinone_pmass1.png' % Nintervals
plotname_two_tiff = endlocation_out+'D_vs_sigma_bulk_better_rms_Nestimates%i_twoinone_pmass1.tiff' % Nintervals
plotname_fit      = endlocation_out+'D_vs_sigma_fit_better_rms_Nestimates%i_pmass1.png' % Nintervals
indfilename       = endlocation_out+'D_vs_sigma_fitindices_better_rms_Nestimates%i_pmass1.txt' % Nintervals

outfile = open(outfilename, 'w')
outfile.write('d   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    psigma = psigmas[i]
    endlocation_in   = 'Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/' % damp +'Pure_bulk/'+ 'Sigma_bead_' +str(psigma) + '/'
    infilename = endlocation_in+'diffusion'+filestext+'_better_rms_Nestimates%i_pmass1.txt' % Nintervals
    metaname   = endlocation_in+'indices_for_fit.txt'
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    infile = open(infilename, "r")
    lines = infile.readlines()
    # Getting the number of lines, etc.
    line = lines[1]
    words = line.split()
    
    # Ds
    DRs[i] = float(words[0])
    Dzs[i] = float(words[2])
    Dparallel[i] = float(words[4])
    # Ds, stdv
    DRs_stdv[i] = float(words[1])
    Dzs_stdv[i] = float(words[3])
    Dparallel_stdv[i] = float(words[5])
    
    infile.close()
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (psigma/2., DRs[i], DRs_stdv[i],  Dzs[i], Dzs_stdv[i], Dparallel[i], Dparallel_stdv[i]))
    
    metafile   = open(metaname, 'r')
    mlines      = metafile.readlines()
    startindex_R  = int(mlines[0].split()[1])
    endindex_R    = int(mlines[1].split()[1])
    startindex_z  = int(mlines[2].split()[1])
    endindex_z    = int(mlines[3].split()[1])
    startindex_p  = int(mlines[4].split()[1])
    endindex_p    = int(mlines[5].split()[1])
    metafile.close()
    indexfile.write('%.2f %i %i %i %i %i %i\n' % (spacing, startindex_R, endindex_R, startindex_z, endindex_z, startindex_p, endindex_p))

outfile.close()

psigmas = np.array(psigmas)
psigmas = psigmas/2.

theory_graph = theory(psigmas)

fig, ax1 = plt.subplots(dpi=300)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.fill_between(psigmas, DRs+DRs_stdv, DRs-DRs_stdv, color='r', alpha=0.3)
ax1.fill_between(psigmas, Dzs+Dzs_stdv, Dzs-Dzs_stdv, color='b', alpha=0.3)
ax1.fill_between(psigmas, Dparallel+Dparallel_stdv, Dparallel-Dparallel_stdv, color='g', alpha=0.3)
ax1.plot(psigmas, DRs, '-s', color='r', label=r'$D_R$')
ax1.plot(psigmas, Dzs, '-o', color='b', label=r'$D_\perp$')
ax1.plot(psigmas, Dparallel, '-*', color='g', label=r'$D_\parallel$')
ax1.plot(psigmas, theory_graph, '--', label=r'$kT/\gamma$')
ax1.set_xlabel(xlabel=r'$a$ (nm)',labelpad=2,fontsize=12)
ax1.set_ylabel(ylabel=r'$D$ (m$^2$/s)',fontsize=12)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
ax1.yaxis.get_offset_text().set_fontsize(12)
ax1.legend(loc='upper left', bbox_to_anchor=(0.07,0.999),fontsize=11)

plt.tight_layout()

ax2 = fig.add_axes([0.4,0.4,0.5,0.5])
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
ax2.loglog(psigmas, DRs, '-s', color='r')
ax2.loglog(psigmas, Dzs, '-o', color='b')
ax2.loglog(psigmas, Dparallel, '-*', color='g')
ax2.loglog(psigmas, theory_graph, '--', label=r'$kT/\gamma$')
ax2.set_xlabel(xlabel=r'$a$ (nm)',labelpad=-1.5,fontsize=12)
ax2.set_ylabel(ylabel=r'$D$ (m$^2$/s)',labelpad=-3,fontsize=12)
ax2.legend(loc='lower left',fontsize=11)
if plot_in_tiff==True:
    plt.savefig(plotname_two_tiff)
else:
    plt.savefig(plotname_two)

plt.show()

## Find exponent:
DRs_all = np.zeros((N,Nintervals))

for i in range(N):
    psigma = psigmas_base[i]
    inlocation = 'Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Pure_bulk/Sigma_bead_' % damp+str(psigma) + '/'
    allDRsname     = inlocation+'diffusion'+filestext+'_better_rms_allDRs_Nestimates%i.txt' % Nintervals
    allDRsfile     = open(allDRsname,'r')
    lines          = allDRsfile.readlines()
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            for j in range(Nintervals):
                DRs_all[i,j] = float(words[j])
    allDRsfile.close()

exponents = np.zeros(Nintervals)
for i in range(Nintervals):
    logds  = np.log(psigmas)
    logDRs = np.log(DRs_all[:,i])
    
    exponents[i] = (logDRs[-1]-logDRs[1])/(logds[-1]-logds[1])
                
exponent_avg, exponent_rms = avg_and_rms(exponents)

print('exponent_avg:',exponent_avg)
print('exponent_rms:',exponent_rms)
