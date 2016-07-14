"""
    Injecting a KV beam into the idealized solenoid field.
"""
from __future__ import division
from warp_init_tools import *   #Not really used butbrings in ParticleDiagnostic.
import numpy as np
import os
import random
from scipy.optimize import newton
from sys import exit
import matplotlib.pyplot as plt



diagDir = 'diags/xySlice/hdf5'

def cleanupPrevious(outputDirectory = diagDir):
    if os.path.exists(outputDirectory):
        files = os.listdir(outputDirectory)
        for file in files:
            if file.endswith('.h5'):
                os.remove(os.path.join(outputDirectory,file))

cleanupPrevious()

#########################
### Initialize Plots  ###
#########################
def setup():
	pass

##########################################
### Create Beam and Set its Parameters ###
##########################################

SC = False # Controls field solve
ptcl_per_step = 20 # number of particles to inject on each step
beam_beta = 0.56823 # v/c for 110 keV electrons


top.lrelativ = True
top.relativity = 1



beam = Species(type=Electron, name='Electron')

if SC == True:
    beam.ibeam = 8e-3 # Set small during solenoid testing


def createKV(npart = ptcl_per_step, vc = beam_beta):
    ptcls = []


    #Set beam size and emittance
    a = 0.05
    b = 0.05 #/ np.sqrt(200.) 
    rmsemitx = 1.e-6 #* 200.
    rmsemity = 1.e-6 
    emitx = 4. * rmsemitx #geometric,total emittance
    emity = 4. * rmsemity
    
    
    #twiss alpha - Currently must be zero to give correct distr.
    alphax = 0
    alphay = 0
    betax = a**2 / emitx
    betay = b**2 / emity
    gammax = (1 + alphax**2) / betax
    gammay = (1 + alphay**2) / betay

    i = 0
    while i < npart:
        x = y = 20 * max([a,b])
        while (x / a)**2 + (y / b)**2 > 1: # a little inefficient but averages ~1.5 cycles 
            x = (1.0 - 2.0 * random.random()) * a
            y = (1.0 - 2.0 * random.random()) * b
        
        R = 1 - (x / a)**2 - (y / b)**2 

        theta = random.random() * 2 * np.pi

        xp = np.cos(theta) * (np.sqrt(R) * emitx) / a
        yp = np.sin(theta) * (np.sqrt(R) * emity) / b

        ptcls.append([x,xp,y,yp,0.0,vc*3e8])
        i += 1
        #print (x / a)**2 + (y / b)**2 + (xp * a / emit)**2 + (yp * b / emit)**2

    return np.array(ptcls)

ptclArray = createKV()

print ptclArray[:,1]
print np.sqrt((ptclArray[:,1] * ptclArray[:,5])**2 + (ptclArray[:,3] * ptclArray[:,5])**2)


def FRBT(array, beta=12.5 /4., alpha=0.0):
    """ 
        Transforms a matched flat beam to a round 'magnetized' beam.
    """

    gamma = (1. - alpha**2) / beta

    R = np.zeros([6,6],dtype='float64')
    R[0,0] = 1. + alpha
    R[0,1] = beta
    R[0,2] = 1. - alpha
    R[0,3] = -beta

    R[1,0] = -gamma
    R[1,1] = 1. - alpha
    R[1,2] = gamma
    R[1,3] = 1. + alpha

    R[2,0] = 1. - alpha
    R[2,1] = -beta
    R[2,2] = 1. + alpha
    R[2,3] = beta

    R[3,0] = gamma
    R[3,1] = 1. + alpha
    R[3,2] = -gamma
    R[3,3] = 1. - alpha

    R[4,4] = 2.
    R[5,5] = 2.

    R = 0.5 * R

    RB = np.zeros_like(array)

    rows = array.shape[0]

    for i in range(rows):
        for j in range(4):
            val = 0
            for k in range(4):
                val += R[j,k] * array[i,k]
            RB[i,j] = val
    RB[:,4] = array[:,4]
    RB[:,5] = array[:,5]

    return RB

FBeam = ptclArray

print "#######################"
print "%s %s %s" % (FBeam[10,1]*FBeam[10,5], FBeam[10,3]*FBeam[10,5], FBeam[10,5])
print "#######################"

def injectelectrons():

    #Change x and y angles to velocities
    beam.addparticles(x=FBeam[:,0],y=FBeam[:,2],z=FBeam[:,4],
    vx=FBeam[:,1] * FBeam[:,5],vy=FBeam[:,3] * FBeam[:,5],vz=FBeam[:,5])
    #beam.addparticles(x=0.08,y=0.0,z=0.0,vx=0.01*3e8,vy=0.,vz=beam_beta*3e8)


#Install injector
installuserinjection(injectelectrons)

derivqty() #Sets addition derived parameters (such as beam.vbeam)



##########################
### Injection Controls ###
##########################


#--- Specify injection of the particles
top.inject      = 6                       # 2 means space-charge limited injection, 6 is user specified
top.npinject    = ptcl_per_step           # Approximate number of particles injected each step
top.ainject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything
top.binject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything


# --- Other injection variables - not sure if these are important 
w3d.l_inj_exact = True                    # if true, position and angle of injected particle are 
                                          # computed analytically rather than interpolated
w3d.l_inj_area  = False            # Not sure what this does




################################
### 3D Simulation Parameters ###
################################

#Set cells
w3d.nx = 64
w3d.ny = 64
w3d.nz = 32


w3d.bound0  = dirichlet
w3d.boundnz = dirichlet
w3d.boundxy = neumann

#Set boundaries
w3d.xmmin = -0.16 
w3d.xmmax =  0.16 
w3d.ymmin = -0.16 
w3d.ymmax =  0.16 
w3d.zmmin =  0.0
w3d.zmmax =  1.0

#top.pbound0 = 0
#top.pboundnz = 0
#top.pboundxy = 0            # Absorbing Boundary for particles


dz =  w3d.zmmax - w3d.zmmin

top.dt = dz / (14000 * beam_beta * 3e8)

top.ibpush   = 2            # set type of pusher to  vXB push without tan corrections
                            ## 0:off, 1:fast, 2:accurate

############################
### Particle Diagnostics ###
############################

diagP0 = ParticleDiagnostic( period=1, top=top, w3d=w3d,
        species= { species.name : species for species in listofallspecies },
        comm_world=comm_world, lparallel_output=False, write_dir = diagDir[:-4] )

diagP = ParticleDiagnostic( period=100, top=top, w3d=w3d,
        species= { species.name : species for species in listofallspecies },
        comm_world=comm_world, lparallel_output=False, write_dir = diagDir[:-4] )

installafterstep( diagP0.write )
installafterstep( diagP.write )

####################################
### Install Ideal Solenoid Field ###
####################################

# Uniform B_z field
bz = np.zeros([w3d.nx,w3d.ny,w3d.nz])
bz[:,:,:] = 1.0
z_start = 0.0
z_stop = 1.0


addnewbgrd(z_start,z_stop,xs=w3d.xmmin,dx=(w3d.xmmax-w3d.xmmin),ys=w3d.ymmin,dy=(w3d.ymmax-w3d.ymmin),
    nx=w3d.nx,ny=w3d.ny,nz=w3d.nz,bz=bz)


# Add a full solenoid
# addnewsolenoid(zi=0.01,
#                zf=0.99,
#                ri=0.16,
#                maxbz=-1.)

#################################
### Generate and Run PIC Code ###
#################################

w3d.solvergeom = w3d.XYZgeom

if SC == True:
    solver = MultiGrid3D()
    registersolver(solver)


package("w3d")
generate()



step(1)
uninstallafterstep( diagP0.write ) # You really don't want 1000 diagnostic files

step(200)

uninstalluserinjection(injectelectrons) # Can uninstall injection to stop

step(8800)


