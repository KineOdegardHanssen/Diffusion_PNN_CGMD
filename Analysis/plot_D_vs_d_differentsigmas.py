import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import gridspec
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

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend', fontsize=12)

def theory(a):
    k = 1.38e-23 # J/K
    T = 310      # K
    kT = k*T
    unittime = 2.38e-11  # s, unit time step
    a0 = 5e-10           # a at sigma = 1
    mass_base = 8.1e-25       # kg, unit mass
    density = 3*mass_base/(4*np.pi*a0**3)
    mass = density*4./3.*np.pi*a**3
    damp = 10*unittime   # damp in SI units
    eta = mass/(6*np.pi*damp*a0)
    return kT/(eta*2*a)

thickness = 500e-9 # m

confignrs      = np.arange(1,1001)
filestext      = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
endlocation1   = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_1/Nocut/'
endlocation2   = 'Diffusion_bead_near_grid/D_vs_d/Brush/Varysigma/'
bulklocation1  = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_1/'
plotname       = endlocation2+'D_vs_d_varysigma_DATA.tiff'
plotname2      = endlocation2+'D_vs_d_varysigma_DATA_Ddiv1.tiff'
plotname_small = endlocation2+'D_vs_d_varysigma_small_DATA.tiff'

bulkfilename = bulklocation1 + 'diffusion_bulk'+filestext
bulkfilename = bulkfilename +'_cut.txt'

## Files to read
brushfilename_dyn  = endlocation1 + 'D_vs_d_better_rms_Nestimates10.txt'

## Reading sigma = 1
# Dynamic sims:
brushfile_dyn = open(brushfilename_dyn, 'r')

lines = brushfile_dyn.readlines()
N_dyn = len(lines)-1

# ds
ds1        = np.zeros(N_dyn)
dz1        = np.zeros(N_dyn)
dz1_rms    = np.zeros(N_dyn)
difft1     = np.zeros(N_dyn)
difft1_rms = np.zeros(N_dyn)

for i in range(1,N_dyn+1):
    words = lines[i].split()
    j = i-1

    dthis         = float(words[0])
    dzthis        = float(words[3])
    dzthis_rms    = float(words[4])
    ds1[j]        = dthis
    dz1[j]        = dzthis
    dz1_rms[j]    = dzthis_rms
    difft1[j]     = thickness**2/(2*dzthis)
    difft1_rms[j] = thickness**2/(2*dzthis**2)*dzthis_rms

brushfile_dyn.close()
ds1        = ds1[5:14]
dz1        = dz1[5:14]
dz1_rms    = dz1_rms[5:14]
difft1     = difft1[5:14]
difft1_rms = difft1_rms[5:14]
N1         = len(ds1)

bulkfilename_II = 'Diffusion_bead_near_grid/D_vs_d/Bulk/Varysigmas/D_vs_sigma_better_rms_Nestimates10.txt'
bulkfile_II = open(bulkfilename_II,'r')
hline = bulkfile_II.readline()
lines = bulkfile_II.readlines()

sigma_walker = []
Dzs_bulk_vm = []
for line in lines:
    words = line.split()
    if len(words)>0:
        sigma_walker.append(float(words[0])/2.)
        Dzs_bulk_vm.append(float(words[3]))
bulkfile_II.close()

mylinestyles = [[1,1,1,1],[2,1,5,1],[3,1,3,1],[10,0],[2,2,10,2],[10,1,1,1]]
mycolors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown'] 

sigma_walker = sigma_walker[:6]

fitds = np.linspace(2,10,100)
sigma_chain  = 1.0
sigma_walker = np.array(sigma_walker)
sigmas = (sigma_walker+sigma_chain)/2.
Nsigmas = len(sigmas)

fig, ax = plt.subplots(1, 1, figsize=(6,5),dpi=300)
fig = plt.figure(figsize=(12,5),dpi=300)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

j=0

# Adding a few measurements:
infilename = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_0.25/Nocut/D_vs_d_sigma0p25.txt'
infile = open(infilename,'r')
lines = infile.readlines()
ds0p25         = []
dz0p25         = []
dz0p25_rms     = []
difft0p25      = []
difft0p25_rms  = []
for line in lines:
    words = line.split()
    if len(words)>0:
        dthis = float(words[0])
        dzthis = float(words[1])
        drmsthis = float(words[2])
        ds0p25.append(dthis)
        dz0p25.append(dzthis)
        dz0p25_rms.append(drmsthis)
        difft0p25.append(thickness**2/(2*dzthis))
        difft0p25_rms.append(thickness**2/(2*dzthis**2)*drmsthis)
    else:
        break
infile.close()
ds0p25 = np.array(ds0p25)
dz0p25 = np.array(dz0p25)
dz0p25_rms = np.array(dz0p25_rms)
difft0p25  = np.array(difft0p25)
difft0p25_rms = np.array(difft0p25_rms)
ax.fill_between(ds0p25, dz0p25+dz0p25_rms, dz0p25-dz0p25_rms, facecolor=mycolors[0], alpha=0.2)
ax.plot(ds0p25, dz0p25, '-o', color=mycolors[0],label=r'$a$ = 0.125 nm')
##
infilename = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_0.5/Nocut/D_vs_d_sigma0p5.txt'
infile = open(infilename,'r')
lines = infile.readlines()
ds0p5     = []
dz0p5     = []
dz0p5_rms = []
difft0p5      = []
difft0p5_rms  = []
for line in lines:
    words = line.split()
    if len(words)>0:
        dthis = float(words[0])
        dzthis = float(words[1])
        drmsthis = float(words[2])
        ds0p5.append(dthis)
        dz0p5.append(dzthis)
        dz0p5_rms.append(drmsthis)
        difft0p5.append(thickness**2/(2*dzthis))
        difft0p5_rms.append(thickness**2/(2*dzthis**2)*drmsthis)
    else:
        break
infile.close()
ds0p5 = np.array(ds0p5)
dz0p5 = np.array(dz0p5)
dz0p5_rms = np.array(dz0p5_rms)
difft0p5  = np.array(difft0p5)
difft0p5_rms = np.array(difft0p5_rms)
ax.fill_between(ds0p5, dz0p5+dz0p5_rms, dz0p5-dz0p5_rms, facecolor=mycolors[1], alpha=0.2)
ax.plot(ds0p5, dz0p5, '-o', color=mycolors[1],label=r'$a$ = 0.25 nm')
##
ax.fill_between(ds1, dz1+dz1_rms, dz1-dz1_rms, facecolor=mycolors[2], alpha=0.2)
ax.plot(ds1, dz1, '-o', color=mycolors[2],label=r'$a$ = 0.5 nm')
##
infilename = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_1.5/Nocut/D_vs_d_sigma1p5.txt'
infile = open(infilename,'r')
lines = infile.readlines()
ds1p5     = []
dz1p5     = []
dz1p5_rms = []
difft1p5      = []
difft1p5_rms  = []
for line in lines:
    words = line.split()
    if len(words)>0:
        dthis = float(words[0])
        dzthis = float(words[1])
        drmsthis = float(words[2])
        ds1p5.append(dthis)
        dz1p5.append(dzthis)
        dz1p5_rms.append(drmsthis)
        difft1p5.append(thickness**2/(2*dzthis))
        difft1p5_rms.append(thickness**2/(2*dzthis**2)*drmsthis)
    else:
        break
infile.close()
ds1p5 = np.array(ds1p5)
dz1p5 = np.array(dz1p5)
dz1p5_rms = np.array(dz1p5_rms)
difft1p5  = np.array(difft1p5)
difft1p5_rms = np.array(difft1p5_rms)
ax.fill_between(ds1p5, dz1p5+dz1p5_rms, dz1p5-dz1p5_rms, facecolor=mycolors[3], alpha=0.2)
ax.plot(ds1p5, dz1p5, '-o', color=mycolors[3],label=r'$a$ = 0.75 nm')
##
infilename = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_2/Nocut/D_vs_d_sigma2.txt'
infile = open(infilename,'r')
lines = infile.readlines()
ds2     = []
dz2     = []
dz2_rms = []
difft2      = []
difft2_rms  = []
for line in lines:
    words = line.split()
    if len(words)>0:
        dthis = float(words[0])
        dzthis = float(words[1])
        drmsthis = float(words[2])
        ds2.append(dthis)
        dz2.append(dzthis)
        dz2_rms.append(drmsthis)
        difft2.append(thickness**2/(2*dzthis))
        difft2_rms.append(thickness**2/(2*dzthis**2)*drmsthis)
    else:
        break
infile.close()
ds2 = np.array(ds2)
dz2 = np.array(dz2)
dz2_rms = np.array(dz2_rms)
difft2  = np.array(difft2)
difft2_rms = np.array(difft2_rms)
ax.fill_between(ds2, dz2+dz2_rms, dz2-dz2_rms, facecolor=mycolors[4], alpha=0.2)
ax.plot(ds2, dz2, '-o', color=mycolors[4],label=r'$a$ = 1 nm')
infilename = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_3/Nocut/D_vs_d_sigma3.txt'
infile = open(infilename,'r')
lines = infile.readlines()
ds3     = []
dz3     = []
dz3_rms = []
difft3      = []
difft3_rms  = []
for line in lines:
    words = line.split()
    if len(words)>0:
        dthis = float(words[0])
        dzthis = float(words[1])
        drmsthis = float(words[2])
        ds3.append(dthis)
        dz3.append(dzthis)
        dz3_rms.append(drmsthis)
        difft3.append(thickness**2/(2*dzthis))
        difft3_rms.append(thickness**2/(2*dzthis**2)*drmsthis)
    else:
        break
infile.close()
ds3 = np.array(ds3)
dz3 = np.array(dz3)
dz3_rms = np.array(dz3_rms)
difft3  = np.array(difft3)
difft3_rms = np.array(difft3_rms)
ax.fill_between(ds3, dz3+dz3_rms, dz3-dz3_rms, facecolor=mycolors[5], alpha=0.2)
ax.plot(ds3, dz3, '-o', color=mycolors[5],label=r'$a$ = 1.5 nm')
ax.set_xlabel('d (nm)', fontsize=12)
ax.set_ylabel(r'$D_\perp$ (m$^2$/s)', fontsize=12)
ax.yaxis.get_offset_text().set_fontsize(12)
plt.legend(loc="upper left",fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tight_layout()
plt.savefig(plotname_small)

xiloc     = 'Diffusion_bead_near_grid/'
xifilename = xiloc+'corrlen_vs_d.txt'

xifile = open(xifilename,'r')
lines = xifile.readlines()

d_xis   = []
xis     = []
for line in lines:
    words = line.split()
    
    d_xi = float(words[0])
    for i in range(N1):
        if d_xi==ds1[i]:
            d_xis.append(d_xi)
            xis.append(float(words[1]))

T   = 310            # K
kB  = 1.380649e-23   # J K-1
eta = 3.611078540740483e-07
factor = kB*T/eta
th0p25 = theory(0.125e-9)
th0p5  = theory(0.25e-9)
th1    = theory(0.5e-9)
th1p5  = theory(0.75e-9)
th2    = theory(1.0e-9)
th3    = theory(1.5e-9)

d0p25_th2 = []
d0p5_th2  = []
d1p5_th2  = []
d2_th2    = []
d3_th2    = []
d0p25_theory = []
d0p5_theory  = []
d1_theory    = []
d1p5_theory  = []
d2_theory    = []
d3_theory    = []
for i in range(N1):
    ## Scaling from Cai
    d0p25_th2.append(dz1[i]/0.25)
    d0p5_th2.append(dz1[i]/0.5)
    d1p5_th2.append(dz1[i]/1.5)
    d2_th2.append(dz1[i]/2.0)
    d3_th2.append(dz1[i]/3.0)
    ## Theory from Cai
    d0p25_theory.append(th0p25)
    d0p5_theory.append(th0p5)
    d1_theory.append(th1)
    d1p5_theory.append(th1p5)
    d2_theory.append(th2)
    d3_theory.append(th3)

j = 0
fitds_stat_store = []
diffusiontimes_stat = []
fig = plt.figure(figsize=(12,5),dpi=300)

gs = gridspec.GridSpec(1, 5)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:5])

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#
ax1.fill_between(ds0p25, dz0p25+dz0p25_rms, dz0p25-dz0p25_rms, facecolor='tab:blue', alpha=0.2)
l1 = ax1.plot(ds0p25, dz0p25,'-o', color='tab:blue',label=r'$a$ = 0.125 nm')
l2 = ax1.plot(ds1, dz1/0.25, '--', color='tab:blue',label=r'$a$ = 0.125 nm, scaled')
#
ax1.fill_between(ds0p5, dz0p5+dz0p5_rms, dz0p5-dz0p5_rms, facecolor='tab:orange', alpha=0.2)
l3 = ax1.plot(ds0p5, dz0p5,'-d', color='tab:orange',label=r'$a$ = 0.25 nm')
l4 = ax1.plot(ds1, dz1/0.5, '--', color='tab:orange',label=r'$a$ = 0.25 nm, scaled')
#
ax1.fill_between(ds1, dz1+dz1_rms, dz1-dz1_rms, facecolor='tab:green', alpha=0.2)
ax1.plot(ds1, dz1,'-^', color='tab:green',label=r'$a$ = 0.5 nm')
#
ax1.fill_between(ds1p5, dz1p5+dz1p5_rms, dz1p5-dz1p5_rms, facecolor='tab:red', alpha=0.2)
l5 = ax1.plot(ds1p5, dz1p5,'-*', color='tab:red',label=r'$a$ = 0.75 nm')
l6 = ax1.plot(ds1, dz1/1.5, '--', color='tab:red',label=r'$a$ = 0.75 nm, scaled')
#
ax1.fill_between(ds2, dz2+dz2_rms, dz2-dz2_rms, facecolor='tab:purple', alpha=0.2)
l7 = ax1.plot(ds2, dz2,'-s', color='tab:purple',label=r'$a$ = 1 nm')
l8 = ax1.plot(ds1, dz1/2., '--', color='tab:purple',label=r'$a$ = 1 nm, scaled')
#
ax1.fill_between(ds3, dz3+dz3_rms, dz3-dz3_rms, facecolor='tab:brown', alpha=0.2)
l9 = ax1.plot(ds3, dz3,'-X', color='tab:brown',label=r'$a$ = 1.5 nm')
l10 = ax1.plot(ds1, dz1/3., '--', color='tab:brown',label=r'$a$ = 1.5 nm, scaled')
#
ax1.set_title('A',loc='left',x=-0.03, y=1.05,fontsize=12)
ax2.set_title('B',loc='left',x=-0.03, y=1.05,fontsize=12)
ax1.set_xlabel(r'$d$ (nm)', fontsize=12)
ax1.set_ylabel(r'$D_\perp$ (m$^2$/s)', fontsize=12)
ax1.yaxis.get_offset_text().set_fontsize(12)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#
ax2.fill_between(ds0p25, difft0p25+difft0p25_rms, difft0p25-difft0p25_rms, facecolor='tab:blue', alpha=0.2)
ax2.plot(ds0p25, difft0p25,'-o', color='tab:blue',label=r'$a$ = 0.125 nm')
#
ax2.fill_between(ds0p5, difft0p5+difft0p5_rms, difft0p5-difft0p5_rms, facecolor='tab:orange', alpha=0.2)
ax2.plot(ds0p5, difft0p5,'-d', color='tab:orange',label=r'$a$ = 0.25 nm')
#
ax2.fill_between(ds1, difft1+difft1_rms, difft1-difft1_rms, facecolor='tab:green', alpha=0.2)
ax2.plot(ds1, difft1,'-^', color='tab:green',label=r'$a$ = 0.5 nm')
#
ax2.fill_between(ds1p5, difft1p5+difft1p5_rms, difft1p5-difft1p5_rms, facecolor='tab:red', alpha=0.2)
ax2.plot(ds1p5, difft1p5,'-*', color='tab:red',label=r'$a$ = 0.75 nm')
#
ax2.fill_between(ds2, difft2+difft2_rms, difft2-difft2_rms, facecolor='tab:purple', alpha=0.2)
ax2.plot(ds2, difft2,'-s', color='tab:purple',label=r'$a$ = 1 nm')
#
ax2.fill_between(ds3, difft3+difft3_rms, difft3-difft3_rms, facecolor='tab:brown', alpha=0.2)
ax2.plot(ds3, difft3,'-X', color='tab:brown',label=r'$a$ = 1.5 nm')
#

ax2.set_xlabel(r'$d$ (nm)', fontsize=12)
ax2.set_ylabel('Diffusion time (s)', fontsize=12)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.yaxis.get_offset_text().set_fontsize(12)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#
ax3.set_frame_on(False)
ax3.get_xaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)
ax3.plot('-o', color='tab:blue',label=r'$a$ = 0.125 nm')
ax3.plot( '--', color='tab:blue',label=r'$a$ = 0.125 nm, scaled')
#
ax3.plot('-d', color='tab:orange',label=r'$a$ = 0.25 nm')
ax3.plot('--', color='tab:orange',label=r'$a$ = 0.25 nm, scaled')
#
ax3.plot('-^', color='tab:green',label=r'$a$ = 0.5 nm')
#
ax3.plot('-*', color='tab:red',label=r'$a$ = 0.75 nm')
ax3.plot('--', color='tab:red',label=r'$a$ = 0.75 nm, scaled')
#
ax3.plot('-s', color='tab:purple',label=r'$a$ = 1 nm')
ax3.plot('--', color='tab:purple',label=r'$a$ = 1 nm, scaled')
#
ax3.plot('-X', color='tab:brown',label=r'$a$ = 1.5 nm')
ax3.plot('--', color='tab:brown',label=r'$a$ = 1.5 nm, scaled')
ax3.legend(loc='lower right')

fig.tight_layout()
plt.savefig(plotname)
plt.show()

############ Divide by Dz ##############################
d1_indices     = np.where(np.isin(ds1,ds0p25))
dz1_small      = dz1[d1_indices]
dz1_rms_small  = dz1_rms[d1_indices]
N0p25 = len(ds0p25)
N0p5 = len(ds0p5)
N1p5 = len(ds1p5)
N2   = len(ds2)
N3   = len(ds3)

dz0p25d1_rms   = np.zeros(N0p25)
dz0p5d1_rms    = np.zeros(N0p5)
dz1p5d1_rms    = np.zeros(N1p5)
dz2d1_rms      = np.zeros(N2)
dz3d1_rms      = np.zeros(N3)
dz0p25d1       = np.zeros(N0p25)
dz0p5d1        = np.zeros(N0p5)
dz1p5d1        = np.zeros(N1p5)
dz2d1          = np.zeros(N2)
dz3d1          = np.zeros(N3)
dz1d0p25       = np.zeros(N0p25)
dz0p5d0p25     = np.zeros(N0p5)
dz1p5d0p25     = []
dz2d0p25       = []
dz3d0p25       = []
dz1d0p25_rms   = np.zeros(N0p25)
dz0p5d0p25_rms = np.zeros(N0p5)
dz1p5d0p25_rms = []
dz2d0p25_rms   = []
dz3d0p25_rms   = []
ds1p5d0p25 = []
ds2d0p25 = []
ds3d0p25 = []
for i in range(N0p25):
    dthis = ds0p25[i]
    for j in range(N1):
        if ds1[j]==dthis:
            dz0p25d1[i]     = dz0p25[i]/dz1[j]
            dz0p5d1[i]      = dz0p5[i]/dz1[j]
            dz1d0p25[i]     = dz1[j]/dz0p25[i]
            dz0p25d1_rms[i] = dz0p25[i]/dz1[j]*np.sqrt((dz0p25_rms[i]/dz0p25[i])**2+(dz1_rms[j]/dz1[j])**2)
            dz0p5d1_rms[i]  = dz0p5[i]/dz1[j]*np.sqrt((dz0p5_rms[i]/dz0p5[i])**2+(dz1_rms[j]/dz1[j])**2)
            dz1d0p25_rms[i] = dz1[j]/dz0p25[i]*np.sqrt((dz0p25_rms[i]/dz0p25[i])**2+(dz1_rms[j]/dz1[j])**2)
    dz0p5d0p25[i]     = dz0p5[i]/dz0p25[i]
    dz0p5d0p25_rms[i] = dz0p5[i]/dz0p25[i]*np.sqrt((dz0p5_rms[i]/dz0p5[i])**2+(dz0p25_rms[i]/dz0p25[i])**2)
for i in range(N1p5):
    dthis = ds1p5[i]
    # Divide by D for sigma=1
    for j in range(N1):
        if ds1[j]==dthis:
            dz1p5d1[i]     = dz1p5[i]/dz1[j]
            dz2d1[i]       = dz2[i]/dz1[j]
            dz1p5d1_rms[i] = dz1p5[i]/dz1[j]*np.sqrt((dz1p5_rms[i]/dz1p5[i])**2+(dz1_rms[j]/dz1[j])**2)
            dz2d1_rms[i]   = dz2[i]/dz1[j]*np.sqrt((dz2_rms[i]/dz2[i])**2+(dz1_rms[j]/dz1[j])**2)
    # Divide by D for sigma=0.25
    for j in range(N0p25):
        if ds0p25[j]==dthis:
            dz1p5d0p25.append(dz1p5[i]/dz0p25[j])
            dz2d0p25.append(dz2[i]/dz0p25[j])
            dz1p5d0p25_rms.append(dz1p5[i]/dz0p25[j]*np.sqrt((dz1p5_rms[i]/dz1p5[i])**2+(dz0p25_rms[j]/dz0p25[j])**2))
            dz2d0p25_rms.append(dz2[i]/dz0p25[j]*np.sqrt((dz2_rms[i]/dz2[i])**2+(dz0p25_rms[j]/dz0p25[j])**2))
            ds1p5d0p25.append(dthis)
            ds2d0p25.append(dthis)
for i in range(N3):
    dthis = ds3[i]
    for j in range(N1):
        if ds1[j]==dthis:
            # Divide by D for sigma=1
            dz3d1[i]     = dz3[i]/dz1[j]
            dz3d1_rms[i] = dz3[i]/dz1[j]*np.sqrt((dz2_rms[i]/dz2[i])**2+(dz1_rms[j]/dz1[j])**2)
    # Divide by D for sigma=0.25
    for j in range(N0p25):
        if ds0p25[j]==dthis:
            dz3d0p25.append(dz3[i]/dz0p25[j])
            dz3d0p25_rms.append(dz3[i]/dz0p25[j]*np.sqrt((dz2_rms[i]/dz2[i])**2+(dz0p25_rms[j]/dz0p25[j])**2))
            ds3d0p25.append(dthis)

dz1p5d0p25     = np.array(dz1p5d0p25)
dz2d0p25       = np.array(dz2d0p25)
dz3d0p25       = np.array(dz3d0p25)
dz1p5d0p25_rms = np.array(dz1p5d0p25_rms)
dz2d0p25_rms   = np.array(dz2d0p25_rms)
dz3d0p25_rms   = np.array(dz3d0p25_rms)
ds1p5d0p25     = np.array(ds1p5d0p25)
ds2d0p25       = np.array(ds2d0p25)
ds3d0p25       = np.array(ds3d0p25)

fig = plt.figure(figsize=(12,5),dpi=300)

gs = gridspec.GridSpec(1, 5)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:5])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#
ax1.fill_between(ds0p25, dz0p25d1+dz0p25d1_rms, dz0p25d1-dz0p25d1_rms, facecolor='tab:blue', alpha=0.2)
ax1.plot(ds0p25, dz0p25d1,'-o', color='tab:blue',label=r'$a$ = 0.125 nm')
ax1.plot(ds0p25,np.zeros(N0p25)+4,'--',color='tab:blue',label=r'$a$ = 0.125 nm, theory')
#
ax1.fill_between(ds0p5, dz0p5d1+dz0p5d1_rms, dz0p5d1-dz0p5d1_rms, facecolor='tab:orange', alpha=0.2)
ax1.plot(ds0p5, dz0p5d1,'-d', color='tab:orange',label=r'$a$ = 0.25 nm')
ax1.plot(ds0p5,np.zeros(N0p25)+2,'--',color='tab:orange',label=r'$a$ = 0.25 nm, theory')
#
ax1.plot(ds0p25, dz1_small/dz1_small,'-^', color='tab:green',label=r'$a$ = 0.5 nm')
#
ax1.fill_between(ds1p5, dz1p5d1+dz1p5d1_rms, dz1p5d1-dz1p5d1_rms, facecolor='tab:red', alpha=0.2)
ax1.plot(ds1p5, dz1p5d1,'-*', color='tab:red',label=r'$a$ = 0.75 nm')
ax1.plot(ds1p5,np.zeros(N0p25)+2./3.,'--',color='tab:red',label=r'$a$ = 0.75 nm, theory')
#
ax1.fill_between(ds2, dz2d1+dz2d1_rms, dz2d1-dz2d1_rms, facecolor='tab:purple', alpha=0.2)
ax1.plot(ds2, dz2d1,'-s', color='tab:purple',label=r'$a$ = 1 nm')
ax1.plot(ds2,np.zeros(N0p25)+0.5,'--',color='tab:purple',label=r'$a$ = 1 nm, theory')
#
ax1.fill_between(ds3, dz3d1+dz3d1_rms, dz3d1-dz3d1_rms, facecolor='tab:brown', alpha=0.2)
ax1.plot(ds3, dz3d1,'-X', color='tab:brown',label=r'$a$ = 1.5 nm')
ax1.plot(ds3,np.zeros(N0p25)+1./3.,'--',color='tab:brown',label=r'$a$ = 1.5 nm, theory')
#
ax1.set_title('A',loc='left',x=-0.03, y=1.05,fontsize=12)
ax2.set_title('B',loc='left',x=-0.03, y=1.05,fontsize=12)
ax1.set_xlabel(r'$d$ (nm)', fontsize=12)
ax1.set_ylabel(r'$D_\perp/D_{\perp(a=0.5\ \mathregular{ nm})}$', fontsize=12)
ax1.yaxis.get_offset_text().set_fontsize(12)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#
ax2.fill_between(ds0p25, difft0p25+difft0p25_rms, difft0p25-difft0p25_rms, facecolor='tab:blue', alpha=0.2)
ax2.plot(ds0p25, difft0p25,'-o', color='tab:blue',label=r'$a$ = 0.125 nm')
#
ax2.fill_between(ds0p5, difft0p5+difft0p5_rms, difft0p5-difft0p5_rms, facecolor='tab:orange', alpha=0.2)
ax2.plot(ds0p5, difft0p5,'-d', color='tab:orange',label=r'$a$ = 0.25 nm')
#
ax2.fill_between(ds1, difft1+difft1_rms, difft1-difft1_rms, facecolor='tab:green', alpha=0.2)
ax2.plot(ds1, difft1,'-^', color='tab:green',label=r'$a$ = 0.5 nm')
#
ax2.fill_between(ds1p5, difft1p5+difft1p5_rms, difft1p5-difft1p5_rms, facecolor='tab:red', alpha=0.2)
ax2.plot(ds1p5, difft1p5,'-*', color='tab:red',label=r'$a$ = 0.75 nm')
#
ax2.fill_between(ds2, difft2+difft2_rms, difft2-difft2_rms, facecolor='tab:purple', alpha=0.2)
ax2.plot(ds2, difft2,'-s', color='tab:purple',label=r'$a$ = 1 nm')
#
ax2.fill_between(ds3, difft3+difft3_rms, difft3-difft3_rms, facecolor='tab:brown', alpha=0.2)
ax2.plot(ds3, difft3,'-X', color='tab:brown',label=r'$a$ = 1.5 nm')
#

ax2.set_xlabel(r'$d$ (nm)', fontsize=12)
ax2.set_ylabel('Diffusion time (s)', fontsize=12)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.yaxis.get_offset_text().set_fontsize(12)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#
ax3.set_frame_on(False)
ax3.get_xaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)
ax3.plot('-o', color='tab:blue',label=r'$a$ = 0.125 nm')
ax3.plot('--', color='tab:blue',label=r'$a$ = 0.125 nm, theory')
#
ax3.plot('-d', color='tab:orange',label=r'$a$ = 0.25 nm')
ax3.plot('--', color='tab:orange',label=r'$a$ = 0.25 nm, theory')
#
ax3.plot('-^', color='tab:green',label=r'$a$ = 0.5 nm')
#
ax3.plot('-*', color='tab:red',label=r'$a$ = 0.75 nm')
ax3.plot('--', color='tab:red',label=r'$a$ = 0.75 nm, theory')
#
ax3.plot('-s', color='tab:purple',label=r'$a$ = 1 nm')
ax3.plot('--', color='tab:purple',label=r'$a$ = 1 nm, theory')
#
ax3.plot('-X', color='tab:brown',label=r'$a$ = 1.5 nm')
ax3.plot('--', color='tab:brown',label=r'$a$ = 1.5 nm, theory')
ax3.legend(loc='lower right')
fig.tight_layout()
plt.savefig(plotname2)
plt.show()
