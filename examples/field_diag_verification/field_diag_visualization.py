"""
    A example of a constant-current beam for verifying the field output
"""
from __future__ import division
import numpy as np
import shutil
from shutil import os
from sys import exit
from warp_init_tools import *   # Not really used but brings in ParticleDiagnostic.
from warp.data_dumping.openpmd_diag import field_diag

from rswarp.diagnostics import FieldDiagnostic
from rswarp.utilities.beam_distributions import createKV

diagDir = 'diags/xySlice/hdf5'
solvertype = []
solvertype += ['magnetostatic']
solvertype += ['electrostatic']
# solvertype += ['electromagnetic']
outputFields = True
# outputFields = False


def cleanupPrevious(outputDirectory=diagDir):
    if os.path.exists(outputDirectory):
        files = os.listdir(outputDirectory)
        for file in files:
            if os.path.isdir(outputDirectory + file):
                shutil.rmtree(outputDirectory + file, ignore_errors=True)
            else:
                os.remove(outputDirectory + file)
                # if file.endswith('.h5'):
                #     os.remove(os.path.join(outputDirectory, file))

cleanupPrevious('diags/')

#########################
### Initialize Plots  ###
#########################

def setup():
    pass

##########################################
### Create Beam and Set its Parameters ###
##########################################

beam_gamma = 116e3/511e3 + 1
beam_beta = np.sqrt(1-1/beam_gamma**2)

top.lrelativ = True
top.relativity = 1

beam = Species(type=Electron, name='e-', fselfb=beam_beta * clight)

beam.ibeam = 10e-6

################################
### 3D Simulation Parameters ###
################################

# Set cells
w3d.nx = 64
w3d.ny = 64
w3d.nz = 32

# Set boundaries
w3d.xmmin = -0.16
w3d.xmmax = 0.16
w3d.ymmin = -0.16
w3d.ymmax = 0.16
w3d.zmmin = 0.0
w3d.zmmax = 1.0

top.pbound0 = absorb
top.pboundnz = absorb
top.pboundxy = absorb

Lz = (w3d.zmmax - w3d.zmmin)
dz =  Lz / w3d.nz
top.dt = (dz) / (beam_beta * clight) / 3  # 3 timesteps to cross a single cell
ptcl_per_step = int(beam.ibeam * top.dt / echarge)  # number of particles to inject on each step

top.ibpush = 2  # 0:off, 1:fast, 2:accurate

# --- Other injection variables
w3d.l_inj_exact = True
w3d.l_inj_area = False

# w3d.solvergeom = w3d.XYZgeom
w3d.solvergeom = w3d.RZgeom

# w3d.bound0 = dirichlet
# w3d.boundnz = dirichlet
w3d.bound0 = periodic
w3d.boundnz = periodic
w3d.boundxy = neumann

package("w3d")
generate()

# fieldperiod = 5e-9 // top.dt  # number of steps between outputting fields
fieldperiod = 20  # number of steps between outputting fields
diags = []
fieldDiags = []
solvers = []

if 'electrostatic' in solvertype:
    if w3d.solvergeom == w3d.XYZgeom:
        solver = MultiGrid3D()
    elif w3d.solvergeom == w3d.RZgeom:
        solver = MultiGrid2D()
    solvers.append(solver)
    fieldDiags += [FieldDiagnostic.ElectrostaticFields(solver=solver, top=top, w3d=w3d, period=fieldperiod)]
if 'magnetostatic' in solvertype:
    solver = MagnetostaticMG()
    solvers.append(solver)
    fieldDiags += [FieldDiagnostic.MagnetostaticFields(solver=solver, top=top, w3d=w3d, period=fieldperiod)]

for solver in solvers:
    solver.mgtol = [0.01] * 3
    registersolver(solver)

if outputFields is True:
    diags += fieldDiags


def generateDist():
    ptclTrans = createKV(
        npart=ptcl_per_step,
        a=0.010,
        b=0.010,
        emitx=4. * 1.e-6,
        emity=4. * 1.e-6
    )

    dv_over_v = 0
    if dv_over_v > 0:
        vzoffset = np.random.normal(scale=dv_over_v, size=ptcl_per_step)
    else:
        vzoffset = np.reshape([0.0] * ptcl_per_step, (ptcl_per_step, 1))

    vz = beam_beta * clight * ([1.0] * ptcl_per_step + vzoffset)
    zemit = dz / 5  # emitting surface 1/5th of a cell forward
    return np.column_stack((ptclTrans, [zemit] * ptcl_per_step, vz))


def injectelectrons():
    ptclArray = generateDist()
    # Change to x and y angles to velocities
    beam.addparticles(
        x=ptclArray[:, 0],
        y=ptclArray[:, 2],
        z=ptclArray[:, 4],
        vx=ptclArray[:, 1] * ptclArray[:, 5],
        vy=ptclArray[:, 3] * ptclArray[:, 5],
        vz=ptclArray[:, 5]
    )

# Install injector
installuserinjection(injectelectrons)

derivqty()  # Sets addition derived parameters (such as beam.vbeam)

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

diagP = ParticleDiagnostic(
    period=1,
    top=top,
    w3d=w3d,
    species={species.name: species for species in listofallspecies},
    comm_world=comm_world,
    lparallel_output=False,
    write_dir=diagDir[:-4])
diags += [diagP]


def writeDiagnostics():
    for d in diags:
        d.write()

installafterstep(writeDiagnostics)

package("w3d")
generate()

step(1)

writeDiagnostics()
diagP.period = 30

step(100)
# stept(0.5e-6)
