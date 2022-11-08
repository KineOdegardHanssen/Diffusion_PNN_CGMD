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

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

def getphi(d,Lzs):
    phi = 1-(4./3*np.pi*(0.5**3)*101./(d**2*Lzs))
    return phi

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
pmass    = 1 
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
endlocation = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'

dynfolder     = 'Diffusion_bead_near_grid/'
statfolder    = 'Diffusion_staticbrush/'
nostifffolder = 'Diffusion_nostiffness/'
if moresigmas==True:
    endlocation = 'Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'
endlocation_nostiff = endlocation + 'Nostiffness/'
endlocation_f = 'Diffusion_forest/D_vs_d/Nocut/'

basepath_base      = 'Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = 'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'
if old_bulk==False:
    bulkfilename  = bulklocation + 'diffusion_bulk'+filestext+'_new.txt'

## Files to read
# D's
brushfilename_dyn      = endlocation + 'D_vs_d_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_f        = endlocation_f + 'D_vs_d_forest_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_stat     = endlocation_static + 'D_vs_d_static_better_rms_Nestimates%i.txt' % Nintervals
brushfilename_nostiff  = endlocation_nostiff + 'D_vs_d_better_rms_Nestimates%i.txt' % Nintervals
# h's
hfilename_dyn       = dynfolder + 'heights_vs_d_dyn.txt'
hfilename_stat      = statfolder + 'heights_vs_d_stat.txt'
hfilename_nostiff   = nostifffolder + 'heights_vs_d_nostiff.txt'

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

hfile_dyn = open(hfilename_dyn,'r')
lines = hfile_dyn.readlines()
Nh_dyn = len(lines)

hspacings_dyn = np.zeros(Nh_dyn)
hs_dyn        = np.zeros(Nh_dyn)
vp_dyn        = np.zeros(Nh_dyn)
for i in range(Nh_dyn):
    words = lines[i].split()
    
    sp = float(words[0])
    h  = float(words[1])
    hspacings_dyn[i] = sp
    hs_dyn[i] = h
    vp_dyn[i] = getphi(sp,h)
    
hfile_dyn.close()

print('spacings_dyn:',spacings_dyn)
print('hspacings_dyn:',hspacings_dyn)

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


hfile_stat = open(hfilename_stat,'r')
lines = hfile_stat.readlines()
Nh_stat = len(lines)

hspacings_stat = np.zeros(Nh_stat)
hs_stat        = np.zeros(Nh_stat)
vp_stat        = np.zeros(Nh_stat)
for i in range(Nh_stat):
    words = lines[i].split()
    
    sp = float(words[0])
    h  = float(words[1])
    hspacings_stat[i] = sp
    hs_stat[i] = h
    vp_stat[i] = getphi(sp,h)
    
hfile_stat.close()

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

vp_f = np.zeros(N_f)
for i in range(1,N_f+1):
    words = lines[i].split()
    j = i-1
    
    sp = float(words[0])
    spacings_f[j] = sp
    DRs_f[j] = float(words[1])
    Dzs_f[j] = float(words[3])
    Dparallel_f[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_f[j] = float(words[2])
    Dzs_stdv_f[j] = float(words[4])
    Dparallel_stdv_f[j] = float(words[6])
    vp_f[j] = getphi(sp,100.0)
    
brushfile_f.close()

######
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

hfile_nostiff = open(hfilename_nostiff,'r')
lines = hfile_nostiff.readlines()
Nh_nostiff = len(lines)

hspacings_nostiff = np.zeros(Nh_nostiff)
hs_nostiff        = np.zeros(Nh_nostiff)
vp_nostiff        = np.zeros(Nh_nostiff)
for i in range(Nh_nostiff):
    words = lines[i].split()
    
    print('words:',words)
    print('h:',h)
    sp = float(words[0])
    h  = float(words[1])
    hspacings_nostiff[i] = sp
    hs_nostiff[i] = h
    vp_nostiff[i] = getphi(sp,h)
    
print('hs_nostiff:',hs_nostiff)
hfile_nostiff.close()

###
#Bulk:

bulkfile  = open(bulkfilename, 'r')
# D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
if old_bulk==True:
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


plt.figure(figsize=(6,5),dpi=300)
plt.fill_between(vp_dyn,Dzs_dyn+Dzs_stdv_dyn,Dzs_dyn-Dzs_stdv_dyn,facecolor='b',alpha=0.2)
plt.fill_between(vp_stat,Dzs_stat+Dzs_stdv_stat,Dzs_stat-Dzs_stdv_stat,facecolor='c',alpha=0.2)
plt.fill_between(vp_f,Dzs_f+Dzs_stdv_f,Dzs_f-Dzs_stdv_f,facecolor='r',alpha=0.2)
plt.fill_between(vp_nostiff,Dzs_nostiff+Dzs_stdv_nostiff,Dzs_nostiff-Dzs_stdv_nostiff,facecolor='lightcoral',alpha=0.2)
plt.plot(vp_dyn,Dzs_dyn,'-o' ,color='b', label=r'Dynamic')
plt.plot(vp_stat,Dzs_stat, '-s',color='c', label=r'Static')
plt.plot(vp_f,Dzs_f, '-v' ,color='r', label=r'Straight')
plt.plot(vp_nostiff,Dzs_nostiff, '-<' ,color='lightcoral', label=r'No stiffness')
plt.legend(loc='upper left')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$D_\perp/D_\mathregular{bulk}$')
plt.savefig('Reviewresponses/D_vs_phi_z.png')


plt.figure(figsize=(6,5),dpi=300)
plt.fill_between(vp_dyn,Dparallel_dyn+Dparallel_stdv_dyn,Dparallel_dyn-Dparallel_stdv_dyn,facecolor='g',alpha=0.2)
plt.fill_between(vp_stat,Dparallel_stat+Dparallel_stdv_stat,Dparallel_stat-Dparallel_stdv_stat,facecolor='limegreen',alpha=0.2)
plt.fill_between(vp_f,Dparallel_f+Dparallel_stdv_f,Dparallel_f-Dparallel_stdv_f,facecolor='black',alpha=0.2)
plt.fill_between(vp_nostiff,Dparallel_nostiff+Dparallel_stdv_nostiff,Dparallel_nostiff-Dparallel_stdv_nostiff,facecolor='lightgray',alpha=0.2)
plt.plot(vp_dyn,Dparallel_dyn, '-*', color='g', label='Dynamic')
plt.plot(vp_stat,Dparallel_stat, '-d', color='limegreen', label='Static')
plt.plot(vp_f,Dparallel_f, '-1', color='black', label='Straight')
plt.plot(vp_nostiff,Dparallel_nostiff, '-X', color='lightgray', label='No stiffness')
plt.legend(loc='upper left')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$D_\parallel/D_\mathregular{bulk}$')
plt.savefig('Reviewresponses/D_vs_phi_parallel.png')
plt.show()

print('phi_dyn:',vp_dyn)
print('phi_stat:',vp_stat)
print('phi_straight:',vp_f)
print('phi_nostiff:',vp_nostiff)

d_dynnostiffdiff = []
dynnostiffdiff = []
dynnostiffdiff_rel = []
d_straightnostiffdiff = []
straightnostiffdiff = []
straightnostiffdiff_rel = []
for i in range(N_dyn):
    print('--------------')
    d = spacings_dyn[i]
    print('d=',d)
    print('phi_dyn:',vp_dyn[i])
    print('phi_stat:',vp_stat[i])
    for j in range(N_f):
        fd = spacings_f[j]
        if fd==d:
            print('phi_straight:',vp_f[j])
    for j in range(N_nostiff):
        nsd = spacings_nostiff[j]
        if nsd==d:
            print('phi_nostiff:',vp_nostiff[j])
            print('DIFF, dyn-nostiff:',vp_dyn[i]-vp_nostiff[j])
            d_dynnostiffdiff.append(d)
            dynnostiffdiff.append(vp_dyn[i]-vp_nostiff[j])
            dynnostiffdiff_rel.append((vp_dyn[i]-vp_nostiff[j])/vp_nostiff[j])
            # Straight to no stiff.
            d_straightnostiffdiff.append(d)
            straightnostiffdiff.append(vp_f[j]-vp_nostiff[j])
            straightnostiffdiff_rel.append((vp_f[j]-vp_nostiff[j])/vp_nostiff[j])

plt.figure(figsize=(6,5))
plt.plot(d_dynnostiffdiff,dynnostiffdiff)
plt.xlabel(r'$d$')
plt.ylabel(r'$\phi_\mathregular{dyn}$-$\phi_\mathregular{nostiff}$')
plt.title(r'Difference in $\phi$ betw. dyn & nostiff systems')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(d_dynnostiffdiff[1:],dynnostiffdiff[1:])
plt.xlabel(r'$d$')
plt.ylabel(r'$\phi_\mathregular{dyn}$-$\phi_\mathregular{nostiff}$')
plt.title(r'Difference in $\phi$ betw. dyn & nostiff systems')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(d_dynnostiffdiff,dynnostiffdiff_rel)
plt.xlabel(r'$d$')
plt.ylabel(r'$(\phi_\mathregular{dyn}$-$\phi_\mathregular{nostiff}$)/$\phi_\mathregular{nostiff}$')
plt.title(r'Relatvive difference in $\phi$ betw. dyn & nostiff systems')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(d_dynnostiffdiff[1:],dynnostiffdiff_rel[1:])
plt.xlabel(r'$d$')
plt.ylabel(r'$(\phi_\mathregular{dyn}$-$\phi_\mathregular{nostiff}$)/$\phi_\mathregular{nostiff}$')
plt.title(r'Relatvive difference in $\phi$ betw. dyn & nostiff systems')
plt.tight_layout()

#### Straight to no stiff

plt.figure(figsize=(6,5))
plt.plot(d_straightnostiffdiff,straightnostiffdiff)
plt.xlabel(r'$d$')
plt.ylabel(r'$\phi_\mathregular{straight}$-$\phi_\mathregular{nostiff}$')
plt.title(r'Difference in $\phi$ betw. straight & nostiff systems')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(d_straightnostiffdiff[1:],straightnostiffdiff[1:])
plt.xlabel(r'$d$')
plt.ylabel(r'$\phi_\mathregular{straight}$-$\phi_\mathregular{nostiff}$')
plt.title(r'Difference in $\phi$ betw. straight & nostiff systems')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(d_straightnostiffdiff,straightnostiffdiff_rel)
plt.xlabel(r'$d$')
plt.ylabel(r'$(\phi_\mathregular{straight}$-$\phi_\mathregular{nostiff}$)/$\phi_\mathregular{nostiff}$')
plt.title(r'Relative difference in $\phi$ betw. straight & nostiff systems')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(d_straightnostiffdiff[1:],straightnostiffdiff_rel[1:])
plt.xlabel(r'$d$')
plt.ylabel(r'$(\phi_\mathregular{straight}$-$\phi_\mathregular{nostiff}$)/$\phi_\mathregular{nostiff}$')
plt.title(r'Relative difference in $\phi$ betw. straight & nostiff systems')
plt.tight_layout()

########

hs_f = np.zeros(N_f)+100.

plt.figure(figsize=(6,5))
plt.plot(spacings_dyn,hs_dyn, '-*', color='g', label='Dynamic')
plt.plot(spacings_stat,hs_stat, '-d', color='limegreen', label='Static')
plt.plot(spacings_f,hs_f, '-1', color='black', label='Straight')
plt.plot(spacings_nostiff,hs_nostiff, '-X', color='lightgray', label='No stiffness')
plt.xlabel(r'$d$')
plt.ylabel(r'$h$')
plt.legend(loc='upper right')
plt.title(r'Spacings vs heights')
plt.tight_layout()

plt.figure(figsize=(6,5))
plt.plot(spacings_dyn,vp_dyn, '-*', color='g', label='Dynamic')
plt.plot(spacings_stat,vp_stat, '-d', color='limegreen', label='Static')
plt.plot(spacings_f,vp_f, '-1', color='black', label='Straight')
plt.plot(spacings_nostiff,vp_nostiff, '-X', color='lightgray', label='No stiffness')
plt.xlabel(r'$d$')
plt.ylabel(r'$\phi$')
plt.legend(loc='lower right')
plt.title(r'Spacings vs porosities')
plt.tight_layout()

plt.show()

