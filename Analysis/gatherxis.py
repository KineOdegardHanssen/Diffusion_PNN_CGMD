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

ds = [1,1.25,1.5,2,3,4,5,6,7,8,10]
N  = len(ds)

sigmas1 = np.zeros(N)
sigmas1p5 = np.zeros(N)
sigmas2 = np.zeros(N)
sigmas3 = np.zeros(N)
for i in range(N):
    sigmas1[i]   = 1
    sigmas1p5[i] = 1.5
    sigmas2[i]   = 2
    sigmas3[i]   = 3

xis_avg = np.zeros(N)
xis_rms = np.zeros(N)

confignrs = np.arange(1,1001)
filestext = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
plotname  = 'd_vs_xi.png'
plotname_wom  = 'd_vs_xi_wom.png'

baseloc     = 'Diffusion_bead_near_grid/'
outfilename = baseloc+'corrlen_vs_d.txt'

outfile     = open(outfilename,'w')

for i in range(N):
    spacing = ds[i]
    print('spacing:',spacing)
    # Dyn
    endlocation = baseloc+'Spacing'+str(spacing)+'/damp10_diffseedLgv/Brush/Sigma_bead_1/Nocut/'
    infilename  = endlocation +'corrlen_'+filestext+'_nocut.txt'
    
    infile = open(infilename,'r')
    line  = infile.readline()
    words = line.split()
    xi_avg = float(words[0])
    xi_rms = float(words[1])
    xis_avg[i] = xi_avg
    xis_rms[i] = xi_rms
    infile.close()
    
    outfile.write('%.2f %.5e %.5e\n' % (spacing,xi_avg,xi_rms))
    print('xi_avg:',xi_avg,'; xi_rms:', xi_rms)
    
outfile.close()

plt.figure(figsize=(6,5),dpi=300)
plt.fill_between(ds, xis_avg+xis_rms, xis_avg-xis_rms, alpha=0.2)
plt.plot(ds,xis_avg,'-o',label=r'$\xi$')
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$\xi$ (nm)')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname)


plt.figure(figsize=(6,5),dpi=300)
plt.fill_between(ds, xis_avg+xis_rms, xis_avg-xis_rms, alpha=0.2)
plt.plot(ds,xis_avg,label=r'$\xi$')
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$\xi$ (nm)')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname_wom)
plt.show()