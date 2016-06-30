"""
    A simple example for defining and loading a user-defined beam distribution and continuously injecting it.
"""
from __future__ import division
from warp_init_tools import *   # Not really used but brings in ParticleDiagnostic.
from warp.particles.ionization import Ionization
from rswarp.diagnostics import FieldDiagnostic
from rswarp.utilities.beam_distributions import createKV
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
                os.remove(os.path.join(outputDirectory, file))

cleanupPrevious()

#########################
### Initialize Plots  ###
#########################
# def setup():
# 	pass
setup(cgmlog=0)

##########################################
### Create Beam and Set its Parameters ###
##########################################

beam_gamma = 116e3/511e3 + 1
beam_beta = np.sqrt(1-1/beam_gamma**2)

top.lrelativ = True
top.relativity = 1

beam = Species(type=Electron, name='e-', fselfb=beam_beta * clight)
h2plus = Species(type=Dihydrogen, charge_state=+1, name='H2+', weight=None)
emittedelec = Species(type=Electron, name='emitted e-', weight=None)

beam.ibeam = 1e-6

################################
### 3D Simulation Parameters ###
################################

#Set cells
w3d.nx = 64
w3d.ny = 64
w3d.nz = 32

#Set boundaries
w3d.xmmin = -0.16
w3d.xmmax = 0.16
w3d.ymmin = -0.16
w3d.ymmax = 0.16
w3d.zmmin = 0.0
w3d.zmmax = 1.0

top.pbound0 = absorb
top.pboundnz = absorb
top.pboundxy = reflect

dz = w3d.zmmax - w3d.zmmin
# top.dt = dz / (1000 * beam_beta * clight)
top.dt = (dz / w3d.nz) / (beam_beta * clight) / 5  # 5 timesteps to cross a single cell
ptcl_per_step = beam.ibeam * top.dt // echarge  # number of particles to inject on each step
print("Timestep is %.2E s" % top.dt)

top.ibpush = 1            # set type of pusher to  vXB push without tan corrections
                            ## 0:off, 1:fast, 2:accurate

# --- Other injection variables - not sure if these are important
w3d.l_inj_exact = True                    # if true, position and angle of injected particle are
w3d.l_inj_area = False            # Not sure what this does

w3d.solvergeom = w3d.XYZgeom

w3d.bound0 = dirichlet
w3d.boundnz = dirichlet
w3d.boundxy = neumann

solver = MagnetostaticMG()
diagF = FieldDiagnostic.MagnetostaticFields(solver, top)
installafterstep( diagF.write )
# solver = MultiGrid3D()
# diagF = FieldDiagnostic.ElectrostaticFields(solver, top)
# installafterstep( diagF.write )
solver.mgtol = [0.01]*3
registersolver(solver)



ptclTrans = createKV(
    npart=ptcl_per_step,
    a=0.010,
    b=0.010,
    emitx=4. * 1.e-6,
    emity=4. * 1.e-6
)

vz = beam_beta * clight
zemit = dz / w3d.nz / 5  # emitting surface 1/5th of a cell forward
ptclArray = np.column_stack((ptclTrans, [zemit] * ptcl_per_step, [vz] * ptcl_per_step))

def injectelectrons():

    #Change to x and y angles to velocities
    beam.addparticles(
        x=ptclArray[:, 0],
        y=ptclArray[:, 2],
        z=ptclArray[:, 4],
        vx=ptclArray[:, 1] * ptclArray[:, 5],
        vy=ptclArray[:, 3] * ptclArray[:, 5],
        vz=ptclArray[:, 5]
    )



#Install injector
installuserinjection(injectelectrons)



####################################
### Ionization of background gas ###
####################################


ioniz = Ionization(stride=100,
    xmin=w3d.xmmin,
    xmax=w3d.xmmax,
    ymin=w3d.ymmin,
    ymax=w3d.ymmax,
    zmin=w3d.zmmin,
    zmax=w3d.zmmax,
    nx=w3d.nx,
    ny=w3d.ny,
    nz=w3d.nz,
    l_verbose=True)

target_density = 1.e20

# ------------ e + H2 -> 2e + H2+


ioniz.add(incident_species=beam,
          emitted_species=[h2plus, emittedelec],
          emitted_energy0=[1, 1], # Array of emission energies (in eV) corresponding to emitted_species
          emitted_energy_sigma=[0, 0.1],
          l_remove_target=False, # Flag for removing target particle
          # Can this be a function of incidence parameters like energy?
          cross_section=4e-23, # Where does this figure come from?
        #   l_verbose=True,
          ndens=target_density)


derivqty() #Sets addition derived parameters (such as beam.vbeam)

##########################
### Injection Controls ###
##########################


# --- Specify injection of the particles
top.inject = 6                       # 2 means space-charge limited injection, 6 is user specified
top.npinject = ptcl_per_step           # Approximate number of particles injected each step
top.ainject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything
top.binject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything







############################
### Particle Diagnostics ###
############################

diagP = ParticleDiagnostic( period=1, top=top, w3d=w3d,
        species= { species.name : species for species in listofallspecies },
        comm_world=comm_world, lparallel_output=False, write_dir = diagDir[:-4] )
installafterstep( diagP.write )


package("w3d")
generate()


step(1)

diagP.write() # You really don't want 1000 diagnostic files
diagP.period = 10

step(49)

# uninstalluserinjection(injectelectrons) # Can uninstall injection to stop

step(250)
