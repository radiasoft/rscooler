"""
    A simple example for defining and loading a user-defined beam distribution and continuously injecting it.
"""
from __future__ import division
from warp_init_tools import *   #Not really used butbrings in ParticleDiagnostic.
import numpy as np
import os
import random
from scipy.optimize import newton
from sys import exit



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

SC = True # Controls field solve
ptcl_per_step = 3000 # number of particles to inject on each step
beam_beta = 0.05 # beta = v / c


top.lrelativ = True
top.relativity = 1



beam = Species(type=Electron, name='Electron')

if SC == True:
    beam.ibeam = 1e-3


def createKV(npart = ptcl_per_step, vc = beam_beta):
    ptcls = []


    #Set beam size and emittance
    a = 0.010
    b = 0.010
    rmsemit = 1.e-6
    emit = 4. * rmsemit #geometric,total emittance
    
    
    #twiss alpha - Currently must be zero to give correct distr.
    alphax = 0
    alphay = 0
    betax = a**2 / emit
    betay = b**2 / emit
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

        xp = np.cos(theta) * (np.sqrt(R) * emit) / a
        yp = np.sin(theta) * (np.sqrt(R) * emit) / b

        ptcls.append([x,xp,y,yp,0.0,vc*3e8])
        i += 1
        #print (x / a)**2 + (y / b)**2 + (xp * a / emit)**2 + (yp * b / emit)**2

    return np.array(ptcls)


ptclArray = createKV()


def injectelectrons():

    #Change to x and y angles to velocities
    beam.addparticles(x=ptclArray[:,0],y=ptclArray[:,2],z=ptclArray[:,4],
    vx=ptclArray[:,1] * ptclArray[:,5],vy=ptclArray[:,3] * ptclArray[:,5],vz=ptclArray[:,5])



#Install injector
installuserinjection(injectelectrons)

derivqty() #Sets addition derived parameters (such as beam.vbeam)



##########################
### Injection Controls ###
##########################


# --- Specify injection of the particles
top.inject      = 6                       # 2 means space-charge limited injection, 6 is user specified
top.npinject    = ptcl_per_step           # Approximate number of particles injected each step
top.ainject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything
top.binject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything


# --- Other injection variables - not sure if these are important 
w3d.l_inj_exact = True                    # if true, position and angle of injected particle are 
                                          # computed analytically rather than interpolated
w3d.l_inj_area         = False            # Not sure what this does




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

top.dt = dz / (1000 * beam_beta * 3e8)

top.ibpush   = 1            # set type of pusher to  vXB push without tan corrections
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

step(49)

uninstalluserinjection(injectelectrons) # Can uninstall injection to stop
step(250)


