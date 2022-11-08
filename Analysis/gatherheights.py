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

def theory(d,A):
    return A*d**(-2./3.)

def fitexp(d,A,b):
    return A*d**(-b)

ds_main  = [1,1.25,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,10,15]
ds_other = [1,2,3,4,5,6,7,8,10,15]

Nmain = len(ds_main)
Nother = len(ds_other)

# Avg
hs_dyn = np.zeros(Nmain)
hs_stat = np.zeros(Nmain)
hs_nostiff = np.zeros(Nother)
# Rms
hs_dyn_rms = np.zeros(Nmain)
hs_stat_rms = np.zeros(Nmain)
hs_nostiff_rms = np.zeros(Nother)

dynloc_base = 'Diffusion_bead_near_grid/'
statloc_base = 'Diffusion_staticbrush/'
nostiffloc_base = 'Diffusion_nostiffness/'

outfilename_dyn     = dynloc_base+'heights_vs_d_dyn.txt'
outfilename_stat    = statloc_base+'heights_vs_d_stat.txt'
outfilename_nostiff = nostiffloc_base+'heights_vs_d_nostiff.txt'

outfile_dyn     = open(outfilename_dyn,'w')
outfile_stat    = open(outfilename_stat,'w')
outfile_nostiff = open(outfilename_nostiff,'w')

for i in range(Nmain):
    spacing = ds_main[i]
    print('spacing:',spacing)
    # Dyn
    endlocation_dyn = dynloc_base+'Spacing'+str(spacing)+'/damp10_diffseedLgv/Brush/Sigma_bead_1/'
    hname_dyn   = endlocation_dyn +'maxz_az_config1to1000.txt'
    
    hfile_dyn = open(hname_dyn,'r')
    line  = hfile_dyn.readline()
    words = line.split()
    h_dyn = float(words[0])
    h_dyn_rms = float(words[1])
    hs_dyn[i] = h_dyn
    hs_dyn_rms[i] = h_dyn_rms
    hfile_dyn.close()
    
    outfile_dyn.write('%.2f %.16e\n' % (spacing,h_dyn))
    
    # Stat
    endlocation_stat = statloc_base+'Spacing'+str(spacing)+'/Radius1/Results/'
    hname_stat   = endlocation_stat +'maxz_az_config1to1000_placements10to11.txt'
    
    hfile_stat = open(hname_stat,'r')
    line  = hfile_stat.readline()
    words = line.split()
    h_stat = float(words[0])
    h_stat_rms = float(words[1])
    hs_stat[i] = h_stat
    hs_stat_rms[i] = h_stat_rms
    hfile_stat.close()
    
    outfile_stat.write('%.2f %.16e\n' % (spacing,h_stat))
outfile_dyn.close()
outfile_stat.close()

for i in range(Nother):
    spacing = ds_other[i]
    print('nostiff, d=',spacing)
    # Dyn
    endlocation_nostiff = nostiffloc_base+'Spacing'+str(spacing)+'/Cut/'
    hname_nostiff   = endlocation_nostiff +'maxz_az_config1to1000.txt'
    
    hfile_nostiff = open(hname_nostiff,'r')
    line  = hfile_nostiff.readline()
    words = line.split()
    h_nostiff = float(words[0])
    h_nostiff_rms = float(words[1])
    hs_nostiff[i] = h_nostiff
    hs_nostiff_rms[i] = h_nostiff_rms
    hfile_nostiff.close()
    
    outfile_nostiff.write('%.2f %.16e\n' % (spacing,h_nostiff))
outfile_nostiff.close()

ds_main = np.array(ds_main)
popt, pcov  = curve_fit(theory, ds_main, hs_dyn, absolute_sigma=True)
Amain      = popt[0]
theorycurve = theory(ds_main,Amain)

ds_other = np.array(ds_other)
popt, pcov  = curve_fit(theory, ds_other, hs_nostiff, absolute_sigma=True)
Aother      = popt[0]
theorycurve_other = theory(ds_other,Aother)

ds_main = np.array(ds_main)
popt, pcov    = curve_fit(theory, ds_main[3:-1], hs_dyn[3:-1], absolute_sigma=True)
Amain_short   = popt[0]
theorycurve_short = theory(ds_main,Amain_short)

ds_main = np.array(ds_main)
popt, pcov    = curve_fit(theory, ds_main[3:-1], hs_stat[3:-1], absolute_sigma=True)
Astat_short   = popt[0]
theorycurve_stat_short = theory(ds_main,Astat_short)

ds_other = np.array(ds_other)
popt, pcov     = curve_fit(theory, ds_other[1:-1], hs_nostiff[1:-1], absolute_sigma=True)
Aother_short   = popt[0]
theorycurve_other_short = theory(ds_other,Aother_short

ds_main = np.array(ds_main)
popt, pcov    = curve_fit(fitexp, ds_main, hs_dyn, absolute_sigma=True)
Amain_short_fe   = popt[0]
bmain_short_fe   = popt[1]
fitexp_short = fitexp(ds_main,Amain_short_fe,bmain_short_fe)


ds_other = np.array(ds_other)
popt, pcov     = curve_fit(fitexp, ds_other[1:-1], hs_nostiff[1:-1], absolute_sigma=True)
Aother_short_fe   = popt[0]
bother_short_fe   = popt[1]
fitexp_other_short = fitexp(ds_other,Aother_short_fe,bother_short_fe)

plt.figure(figsize=(6,5),dpi=300)
plt.plot(ds_main,hs_dyn,label='Dyn.')
plt.plot(ds_main,hs_stat,label='Stat.')
plt.plot(ds_other,hs_nostiff,label='No stiff.')
plt.plot(ds_main,theorycurve_short,'--',color='darkblue',label='Theory (Dyn.)')
plt.plot(ds_main,theorycurve_stat_short,'--',color='darkorange',label='Theory (Stat.)')
plt.plot(ds_other,theorycurve_other_short,'--',color='darkgreen',label='Theory (No stiff.)')
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$h$ (nm)')
plt.title(r'$h$ vs $d$, different brush types')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5),dpi=300)
plt.fill_between(ds_main, hs_dyn+hs_dyn_rms, hs_dyn-hs_dyn_rms, facecolor='tab:blue', alpha=0.2)
plt.plot(ds_main,hs_dyn, color='tab:blue', label='Dyn.')
plt.fill_between(ds_main, hs_stat+hs_stat_rms, hs_stat-hs_stat_rms, facecolor='tab:orange', alpha=0.2)
plt.plot(ds_main,hs_stat, color='tab:orange',label='Stat.')
plt.fill_between(ds_other, hs_nostiff+hs_nostiff_rms, hs_nostiff-hs_nostiff_rms, facecolor='tab:green', alpha=0.2)
plt.plot(ds_other,hs_nostiff, color='tab:green', label='No stiff.')
plt.plot(ds_main,theorycurve_short,'--',color='darkblue',label='Theory (Dyn.)')
plt.plot(ds_main,theorycurve_stat_short,'--',color='chocolate',label='Theory (Stat.)')
plt.plot(ds_other,theorycurve_other_short,'--',color='darkgreen',label='Theory (No stiff.)')
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$L_z$ (nm)')
plt.axis([1,15,25,240])
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()

print('Amain:',Amain_short)
print('Astat:',Astat_short)
print('Aother:',Aother_short)
print('ds_main[:-4]:',ds_main[:-4])
print('bmain_short_fe',bmain_short_fe)
print('bother_short_fe',bother_short_fe)