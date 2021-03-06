import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

###### Tricky parameter: ##################################################
long = False # Only set long=True if you have long runs for the selected d's
###########################################################################

# Input parameters for file selection:
bulkdiffusion = False
substrate     = False
cutit         = False

spacings = [3,3.5,4]
psigma   = 1
damp     = 10

Nsteps = 2001
confignrs = np.arange(1,1001)

# Make file names

for spacing in spacings:
    if bulkdiffusion==True:
        parentfolder = 'Pure_bulk/'
        filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
        if substrate==True:
            parentfolder = 'Bulk_substrate/'
    else:
        parentfolder = 'Brush/'
        filestext     = '_config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_hex'
    filestext = filestext + '_nocut'
    endlocation   = '/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/' % damp +parentfolder+ 'Sigma_bead_' +str(psigma) + '/Hexagonal/Nocut/'
    infilename      = endlocation+'av_ds'+filestext+'.txt'  
    if long==True:
        infilename  = endlocation+'av_ds'+filestext+'_long.txt' 
    outfilename     = endlocation+'diffusion'+filestext+'.txt'
    metaname        = endlocation+'diffusion_metadata'+filestext+'.txt'
    plotname_R      = endlocation+'diffusion_R'+filestext+'.png'
    plotname_dz     = endlocation+'diffusion_dz'+filestext+'.png'
    plotname_dpar   = endlocation+'diffusion_dpar'+filestext+'.png'
    
    # Choosing which part to fit: in file
    rangefilename = endlocation+'indices_for_fit.txt'
    rangefile     = open(rangefilename,'r')
    lines         = rangefile.readlines()
    
    startindex_R   = int(lines[0].split()[1])
    endindex_R     = int(lines[1].split()[1])
    startindex_ort = int(lines[2].split()[1])
    endindex_ort   = int(lines[3].split()[1])
    startindex_par = int(lines[4].split()[1])
    endindex_par   = int(lines[5].split()[1])
    rangefile.close()

    startindices = np.array([startindex_R, startindex_ort,startindex_par])
    endindices   = np.array([endindex_R, endindex_ort, endindex_par])

    startindex = min(startindices)
    endindex   = max(endindices)

    plotseed = 0
    plotdirs = False
    test_sectioned = False
    unitlength   = 1e-9
    unittime     = 2.38e-11 # s
    timestepsize = 0.00045*unittime

    # Set up arrays
    Nsteps_R   = endindex_R-startindex_R
    Nsteps_ort = endindex_ort-startindex_ort 
    Nsteps_par = endindex_par-startindex_par 
    times_R    = np.zeros(Nsteps_R)
    times_ort  = np.zeros(Nsteps_ort)
    times_par  = np.zeros(Nsteps_par)
    dR2s   = np.zeros(Nsteps_R) 
    dz2s   = np.zeros(Nsteps_ort) 
    dpar2s = np.zeros(Nsteps_par)

    # Read file
    infile = open(infilename,'r')
    lines = infile.readlines()

    iR = 0
    iort = 0
    ipar = 0

    #    	0			1		2		3		4			5		6
    # (times_single[i], times_single_real[i], averageRs_SI[i], averagedxs_SI[i], averagedys_SI[i], averagedzs_SI[i], averagedparallel_SI[i]))
    #Loop over sections
    for i in range(startindex,endindex):
        j     = i-startindex
        line  = lines[i]
        words = line.split()
        if i>=startindex_R and i<endindex_R:
            times_R[iR]  = float(words[1])
            dR2s[iR]     = float(words[2])
            iR+=1
        if i>=startindex_ort and i<endindex_ort:
            times_ort[iort] = float(words[1])
            dz2s[iort]      = float(words[5])
            iort+=1
        if i>=startindex_par and i<endindex_par:
            times_par[ipar] = float(words[1])
            dpar2s[ipar]    = float(words[6])
            ipar+=1

    infile.close()

    '''
    plt.figure(figsize=(6,5))
    plt.plot(times_R, dR2s, label='Average, brush')
    #plt.plot(times_R, fit_R2, '--', label='Fit, average, brush')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'Distance$^2$ [in unit length]')
    plt.title('RMSD in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    #plt.savefig(plotname_R)
    plt.show()
    '''

    # line fit, all
    coeffs, covs = polyfit(times_R, dR2s, 1, full=False, cov=True)
    a_R2 = coeffs[0]
    b_R2 = coeffs[1]
    D_R2 = a_R2/6.
    rms_D_R2 = np.sqrt(covs[0,0])/6.
    rms_b_R2 = np.sqrt(covs[1,1])
     
    fit_R2 = a_R2*times_R+b_R2

    # line fit, dz
    coeffs, covs = polyfit(times_ort, dz2s, 1, full=False, cov=True)
    a_z2 = coeffs[0]
    b_z2 = coeffs[1]
    D_z2 = a_z2/2.
    rms_D_z2 = np.sqrt(covs[0,0])/2.
    rms_b_z2 = np.sqrt(covs[1,1])

    fit_z2 = a_z2*times_ort+b_z2

    # line fit, parallel
    coeffs, covs = polyfit(times_par, dpar2s, 1, full=False, cov=True)
    a_par2     = coeffs[0]
    b_par2     = coeffs[1]
    D_par2     = a_par2/4.
    rms_D_par2 = np.sqrt(covs[0,0])/4.
    rms_b_par2 = np.sqrt(covs[1,1])

    fit_par2 = a_par2*times_par+b_par2

    outfile = open(outfilename,'w')
    outfile.write('D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (D_R2, rms_D_R2, b_R2, rms_b_par2, D_z2, rms_D_z2, b_z2, rms_b_par2, D_par2, rms_D_par2, b_par2, rms_b_par2))
    outfile.close()

    metafile = open(metaname, 'w')
    metafile.write('startindex: %i\n' % startindex)
    metafile.write('endindex: %i\n' % endindex)
    metafile.write('startindex_R: %i\n' % startindex_R)
    metafile.write('endindex_R: %i\n' % endindex_R)
    metafile.write('startindex_ort: %i\n' % startindex_ort)
    metafile.write('endindex_ort: %i\n' % endindex_ort)
    metafile.write('startindex_par: %i\n' % startindex_par)
    metafile.write('endindex_par: %i\n' % endindex_par)
    if long==True:
        metafile.write('long file used')
    metafile.close()

    plt.figure(figsize=(6,5))
    plt.plot(times_R, dR2s, label='Average, brush')
    plt.plot(times_R, fit_R2, '--', label='Fit, average, brush')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'Distance$^2$ [in unit length]')
    plt.title('RMSD in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.savefig(plotname_R)

    plt.figure(figsize=(6,5))
    plt.plot(times_ort, dz2s, label='Average, brush')
    plt.plot(times_ort, fit_z2, '--', label='Fit, average, brush')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'$dz^2$ [in unit length]')
    plt.title('$dz^2$ in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.savefig(plotname_dz)

    plt.figure(figsize=(6,5))
    plt.plot(times_par, dpar2s, label='Average, brush')
    plt.plot(times_par, fit_par2, '--', label='Fit, average, brush')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'$dx^2+dy^2$ [in unit length]')
    plt.title('$dx^2+dy^2$ in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.savefig(plotname_dpar)

    print('spacing:', spacing)
    print('psigma:', psigma)
