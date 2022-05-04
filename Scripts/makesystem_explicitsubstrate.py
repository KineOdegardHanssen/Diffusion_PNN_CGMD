import matplotlib.pyplot as plt               # To plot
import numpy as np
import random
import math

onechain = False

# Set variables
N             = 100 # Number of bond vectors = number of units - 1 # Per chain
if onechain==True:
    M             = 1   # Number of chains
    L2            = 1   # Determining the shape of the L1xL2 polymer grid on the surface
    gridspacing   = 2   # This will affect the box size
else:
    M             = 9   # Number of chains
    L2            = 3   # Determining the shape of the L1xL2 polymer grid on the surface
    gridspacing   = 4.5 # The spacing in our grid
substratecharge = 0
charge          = -1
mass            = 1
Nelem           = M*(N+1)
Nbonds          = M*N
Nangles         = M*(N-1)
xpos            = []
ypos            = []
zpos            = []
atomtypes_all   = []
molID           = []
qs              = []
vxs             = []
vys             = []
vzs             = []
vectors         = []
positions       = []

substrate_set = np.zeros((M,2))

for j in range(M):
    # Starting the chain
    x_accm = (j//L2)*gridspacing
    y_accm = j%L2*gridspacing
    z_accm = 0
    substrate_set[j,0] = x_accm
    substrate_set[j,1] = y_accm  
    positions.append([x_accm,y_accm,0])
    xpos.append(x_accm)
    ypos.append(y_accm)
    zpos.append(0)
    vxs.append(np.random.normal(0,1))
    vys.append(np.random.normal(0,1))
    vzs.append(np.random.normal(0,1))
    qs.append(charge)
    molID.append(j+1)
    atomtypes_all.append(1)
    
    # We want another fixed atom per chain:
    chain_twofirst = 0
    # Set up system
    # (Generate random configuration to start with)
    for i in range(N):
        x = 0
        y = 0
        z = 1
        x_accm += x
        y_accm += y
        z_accm += z
        qs.append(charge)
        xpos.append(x_accm)
        ypos.append(y_accm)
        zpos.append(z_accm)
        # All atoms are of type 2 (=moving):
        
        # Except the first two:
        if chain_twofirst==0:
            atomtypes_all.append(1)
            chain_twofirst = 1 # This test is only true once per chain
        else:
            atomtypes_all.append(2)
        molID.append(j+1)
        # Set velocity
        vxs.append(np.random.normal(0,1))
        vys.append(np.random.normal(0,1))
        vzs.append(np.random.normal(0,1))
        qs.append(charge)

### Data on system
xmax = max(xpos)+0.5*gridspacing
ymax = max(ypos)+0.5*gridspacing
zmax = max(zpos)*3.0
zmin = min(zpos)*3.0
xmin = min(xpos)-0.5*gridspacing
ymin = min(ypos)-0.5*gridspacing

# Now I know how big the system is. Next: Setting the substrate.
d    = 1 # One nm (=sigma)
Lx   = xmax-xmin
Ly   = ymax-ymin
Lz   = zmax-zmin
Nx   = int(math.floor(Lx/d))
Ny   = int(math.floor(Ly/d))
breakit = False
Nall = Nelem
for i in range(Nx):
    for j in range(Ny):
        x = i*d+xmin
        y = j*d+ymin
        for k in range(M):
            if substrate_set[k,0]==x and substrate_set[k,1]==y:
                breakit = True     # Do not set a substrate bead if the chain is tethered to this point.
                continue
        if breakit==True:
            breakit = False
            continue
        atomtypes_all.append(4)    # Reserve 3 for free bead
        qs.append(substratecharge) # 
        xpos.append(x)
        ypos.append(y)
        zpos.append(0)
        molID.append(M+2)   # Reserve M+1 for free bead
        vxs.append(np.random.normal(0,1))
        vys.append(np.random.normal(0,1))
        vzs.append(np.random.normal(0,1))
        Nall += 1
        

### Print to file
if mass<1:
    outfilename = 'data.chaingrids_substrate_N%i_Nchains%i_Ly%i_gridspacing'% (Nelem,M,L2)+str(gridspacing)+'_twofixed_charge%i_mass%.5f' % (charge,mass)
else:
    outfilename = 'data.chaingrids_substrate_N%i_Nchains%i_Ly%i_gridspacing'% (Nelem,M,L2)+str(gridspacing)+'_twofixed_charge%i_mass%i' % (charge,mass)
outfile = open(outfilename,'w')
outfile.write('LAMMPS data file via python scripting, version 11 Aug 2017, timestep = 0\n\n%i atoms \n%i bonds \n4 atom types\n1 bond types \n%i angles \n1 angle types\n\n' % (Nall,Nbonds, Nangles))
outfile.write('%.16e %.16e xlo xhi\n' % (xmin, xmax))
outfile.write('%.16e %.16e ylo yhi\n' % (ymin, ymax))
outfile.write('%.16e %.16e zlo zhi\n\n' % (zmin, zmax))
outfile.write('Masses\n\n')
outfile.write('1 %.5f\n2 %.5f\n3 %.5f\n4 %.5f\n' % (mass,mass,mass,mass))

outfile.write('\nAtoms # full\n\n')

for i in range(Nall):
    # For atom_type full:
    # atom-ID molecule-ID atom-type q x y z
    outfile.write('%i %i %i %i %.16e %.16e %.16e\n' % ((i+1), molID[i], atomtypes_all[i], qs[i], xpos[i], ypos[i], zpos[i]))

outfile.write('\nVelocities\n\n')
for i in range(Nall):
    outfile.write('%i %.16e %.16e %.16e\n' % ((i+1), vxs[i], vys[i], vzs[i]))

outfile.write('\nBonds\n\n')

counter = 0
for j in range(M):
    for i in range(N):
        counter += 1
        k = j*(N+1)+i+1
        outfile.write('%i 1 %i %i\n' % (counter, k, (k+1)))

outfile.write('\nAngles\n\n')

counter = 0
for j in range(M):
    for i in range(N-1):
        counter += 1
        k = j*(N+1)+i+1
        outfile.write('%i 1 %i %i %i\n' % (counter, k, (k+1), (k+2)))

outfile.close()

