import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

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

# Input parameters for file selection: 
psigma   = 1 # For instance 
spacings = [1,2,3,4,5,6,7,8,10,15,25,50,75,100]
damp     = 10
N        = len(spacings)
# Input booleans for file selection:
bulkdiffusion = False
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
    parentfolder = 'Pure_bulk/'
    filestext    = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
    systemtype   = 'bulk'
    if substrate==True:
        parentfolder = 'Bulk_substrate/'
        systemtype   = 'substrate'
else:
    parentfolder = 'Brush/'
    systemtype   = 'brush'
    filestext    = '_config'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation_out = '/Diffusion_forest/D_vs_d/Nocut/'
outfilename  = endlocation_out+'D_vs_d_forest.txt'
plotname     = endlocation_out+'D_vs_d_forest.png'
plotname_fit = endlocation_out+'D_vs_d_fit_forest.png'
indfilename  = endlocation_out+'D_vs_d_fitindices_forest.txt'

outfile = open(outfilename, 'w')
outfile.write('d   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    spacing = spacings[i]
    endlocation_in   = '/Diffusion_forest/Spacing'+str(spacing)+'/Results/'
    infilename = endlocation_in+'diffusion'+filestext+'_nocut.txt'
    metaname   = endlocation_in+'diffusion_metadata'+filestext+'_nocut.txt'
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    infile = open(infilename, "r")
    lines = infile.readlines() # This takes some time
    # Getting the number of lines, etc.
    line = lines[1]
    words = line.split()
    
    # Ds
    DRs[i] = float(words[0])
    Dzs[i] = float(words[4])
    Dparallel[i] = float(words[8])
    # Ds, stdv
    DRs_stdv[i] = float(words[1])
    Dzs_stdv[i] = float(words[5])
    Dparallel_stdv[i] = float(words[9])
    
    # bs
    bRs[i] = float(words[2])
    bzs[i] = float(words[6])
    bparallel[i] = float(words[10])
    
    # bs, stdv
    bRs_stdv[i] = float(words[3])
    bzs_stdv[i] = float(words[7])
    bparallel_stdv[i] = float(words[11])
    
    infile.close()
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n' % (spacing, DRs[i], DRs_stdv[i], bRs[i], bRs_stdv[i], Dzs[i], Dzs_stdv[i], bzs[i], bzs_stdv[i], Dparallel[i], Dparallel_stdv[i], bparallel[i], bparallel[i]))
    
    metafile = open(metaname, 'r')
    mlines   = metafile.readlines()
    startindex_R   = int(mlines[0].split()[1])
    endindex_R     = int(mlines[1].split()[1])
    startindex_ort = int(mlines[2].split()[1])
    endindex_ort   = int(mlines[3].split()[1])
    startindex_par = int(mlines[4].split()[1])
    endindex_par   = int(mlines[5].split()[1])
    metafile.close()
    indexfile.write('%i %i %i %i %i %i\n' % (startindex_R, endindex_R, startindex_ort, endindex_ort, startindex_par, endindex_par))

outfile.close()

spacings = np.array(spacings)
    
plt.figure(figsize=(6,5))
plt.errorbar(spacings, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.errorbar(spacings, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.errorbar(spacings, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$d$')
plt.ylabel(r'Diffusion constant $D$')
plt.title('Diffusion constant $D$, d = %i nm, dynamic %s' % (spacing, systemtype))
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname)

