"""
    Testing gyromotion in uniform magnetic field with RZ solver.
"""
from __future__ import division
from warp_init_tools import * 
import numpy as np
import os
import random
from scipy.optimize import newton
from sys import exit
import sys
import pickle
from datetime import datetime
sys.path.append('/Users/chall/research/github/rswarp/')
from rswarp.diagnostics import FieldDiagnostic




diagDir = 'diags/rz/hdf5'
diagFDir = ['rzdiags/rz/fields/magnetic','diags/rz/fields/electric']

def cleanupPrevious(outputDirectory = diagDir, fieldDirectory = diagFDir):
    if os.path.exists(outputDirectory):
        files = os.listdir(outputDirectory)
        for file in files:
            if file.endswith('.h5'):
                os.remove(os.path.join(outputDirectory,file))

    for directory in fieldDirectory:
        if os.path.exists(directory):
            files = os.listdir(directory)
            for file in files:
                if file.endswith('.h5'):
                    os.remove(os.path.join(directory,file))

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

ptcl_per_step = 80000 # number of particles to inject on each step

beam_beta = 0.56823 # v/c for 110 keV electrons
beam_current = 10e-6

beam_weight = 0.5 * beam_current / (echarge * beam_beta * clight) / ptcl_per_step 



top.lrelativ = True
top.relativity = 1

beam = Species(type=Electron, name='Electron', weight=beam_weight)

if SC == False:
    beam.sw = 0.0 # Turn off SC


def generateDist():
    ptclTrans = createKV(
        npart=ptcl_per_step,
        a=0.010,
        b=0.010,
        emitx=4. * 1.e-6,
        emity=4. * 1.e-6
    )

    zrand = np.random.rand(ptcl_per_step,)
    zvel = np.ones_like(zrand) * beam_beta * clight

    return np.column_stack((ptclTrans, zrand, zvel))


def createmybeam():
    ptclArray = generateDist()
    beam.addparticles(x=ptclArray[:,0],y=ptclArray[:,2],z=ptclArray[:,4],
    vx=ptclArray[:,1] * ptclArray[:,5],vy=ptclArray[:,3] * ptclArray[:,5],vz=ptclArray[:,5])



derivqty() #Sets addition derived parameters (such as beam.vbeam)


################################
### 3D Simulation Parameters ###
################################

#Set cells
w3d.nx = 64
w3d.ny = 64
w3d.nz = 32

w3d.bound0 = periodic
w3d.boundnz = periodic

# w3d.bound0  = dirichlet
# w3d.boundnz = dirichlet

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

##########################
### Injection Controls ###
##########################


#--- Specify injection of the particles
top.inject      = 6                       # 2 means space-charge limited injection, 6 is user specified
top.npinject    = total_injected_particles * top.dt * clight * beam_beta / dz # Approximate number of particles injected each step
                                                                              # or average number of particles in interval of a step
                                                                              # will determine current if ibeam is set and beam.sw = 0


top.ainject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything
top.binject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything


# --- Other injection variables - not sure if these are important 
w3d.l_inj_exact = True                    # if true, position and angle of injected particle are 
                                          # computed analytically rather than interpolated
w3d.l_inj_area  = False                   # Not sure what this does


############################
### Particle Diagnostics ###
############################

diagP0 = ParticleDiagnostic( period=1, top=top, w3d=w3d,
        species= { species.name : species for species in listofallspecies },
        comm_world=comm_world, lparallel_output=False, write_dir = diagDir[:-4] )


installafterstep( diagP0.write )


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


#################################
### Generate and Run PIC Code ###
#################################


w3d.solvergeom = w3d.RZgeom

package("w3d") # package/generate Must be called after geometry is set 
generate()     # 


if SC == True:
    solverB = MagnetostaticMG()
    registersolver(solverB)
    solverB.mgtol = [0.01] * 3
    registersolver(solverB)
    solverE = MultiGrid2D()
    registersolver(solverE)


installparticleloader(createmybeam) # for particleloader the call Must be between 1st and 2nd generate calls (or beam doubles)


package("w3d") # package/generate must be called a second time after solver set
generate()     #

fieldperiod = 500
diagBF = FieldDiagnostic.MagnetostaticFields(solver=solverB, top=top, w3d=w3d, period=fieldperiod)
installafterstep(diagBF.write)
diagEF = FieldDiagnostic.ElectrostaticFields(solver=solverE, top=top, w3d=w3d, period=fieldperiod)
installafterstep(diagEF.write)


############
simulation_parameters = {
    'space_charge' : SC,
    'magnetostatic_solver' : SC, # Currently always running solver if SC on
    'electrostatic_solver' : False, # Not being used right now
    'solver_geometry' : w3d.solvergeom,
    'grid_nodes' : (w3d.nx,w3d.ny,w3d.nz),
    'z_boundary_condition' : (w3d.bound0, w3d.boundnz),
    'xy_boundary_condition' : w3d.boundxy,
    'timestep' : top.dt,
    'beam_current' : beam.ibeam,
    'total_particles_injected' : ptclArray.shape[0],
    'run_start_time' : format(datetime.now())
}

pickle.dump(simulation_parameters, open("simulation_parameters.p", 'wb'))
#############


step(1)
diagP0.period = 100 
step(5000)




