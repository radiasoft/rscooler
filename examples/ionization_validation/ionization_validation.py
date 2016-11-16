"""
    A example demonstrating ionization of background H2 gas by an electron beam
"""
from __future__ import division
import numpy as np
import shutil
import time
from shutil import os

import warpoptions
warpoptions.ignoreUnknownArgs = True

import argparse
from warp import *
from warp.data_dumping.openpmd_diag.particle_diag import ParticleDiagnostic
from warp.data_dumping.openpmd_diag import field_diag

from rswarp.diagnostics import FieldDiagnostic
from rswarp.utilities.beam_distributions import createKV
from rswarp.utilities.ionization import Ionization

from warp.diagnostics import gistdummy as gist

import rsoopic.h2crosssections as h2crosssections
from rsoopic.h2crosssections import h2_ioniz_crosssection

solvertype = []
# no field solve
# top.fstype = -1
# solvertype += ['magnetostatic']
solvertype += ['electrostatic']
# outputFields = True
outputFields = False
fieldperiod = 100  # number of steps between outputting fields
particleperiod = 100  # number of steps between outputting particles

solvertolerance = 0.01

simulateIonization = True
# simulateIonization = False

if simulateIonization is True:
    diagDir = 'diags.with/'
else:
    diagDir = 'diags.without/'


def cleanupPrevious(outputDirectory=diagDir):
    if os.path.exists(outputDirectory):
        shutil.rmtree(outputDirectory, ignore_errors=True)

cleanupPrevious()

#####################################################
# Initialize Warp output (mostly disabling it)      #
#####################################################

top.lprntpara = False  # Do not print parameters when generating w3d
top.lpsplots = top.never  # Never generate plots

##########################################
# Create Beam and Set its Parameters     #
##########################################

parser = argparse.ArgumentParser()
parser.add_argument(
    '--beamke',
    type=float,
    action='store',
    dest='beamke',
    help='Electon beam kinetic energy (eV)',
    default=1000
)
parser.add_argument(
    '--dt',
    type=float,
    action='store',
    dest='dt',
    help='Time step (s)',
    default=0.1e-9
)
parser.add_argument(
    '--time',
    type=float,
    action='store',
    dest='time',
    help='Simulation time (s)',
    default=100e-9
)
args = parser.parse_args()
beam_ke = args.beamke  # beam kinetic energy, in eV
beam_gamma = beam_ke/511e3 + 1
beam_beta = np.sqrt(1-1/beam_gamma**2)
diagDir = ("%.3fkeV" % (beam_ke/1e3)) + diagDir

top.lrelativ = True
top.relativity = 1

sw = 100

# beam = Species(type=Electron, name='e-', fselfb=beam_beta * clight, weight=sw)
beam = Species(type=Electron, name='e-', weight=0)
h2plus = Species(type=Dihydrogen, charge_state=+1, name='H2+', weight=0)
emittedelec = Species(type=Electron, name='emitted e-', weight=0)

beam.ibeam = 0.1e-3

################################
# 3D Simulation Parameters     #
################################

# Set cells
w3d.nx = 32
w3d.ny = 32
w3d.nz = 64

# Set boundaries
w3d.xmmin = -0.16
w3d.xmmax = 0.16
w3d.ymmin = -0.16
w3d.ymmax = 0.16
w3d.zmmin = 0.00
w3d.zmmax = 0.20

top.pbound0 = absorb
top.pboundnz = absorb
top.pboundxy = absorb

Lz = (w3d.zmmax - w3d.zmmin)
dz =  Lz / w3d.nz
# top.dt = (dz) / (beam_beta * clight) / 3  # 3 timesteps to cross a single cell
top.dt = args.dt
ptcl_per_step = int(beam.ibeam * top.dt / echarge / sw)  # number of particles to inject on each step

top.ibpush = 1  # 0:off, 1:fast, 2:accurate

# --- Other injection variables
w3d.l_inj_exact = True
w3d.l_inj_area = False

w3d.solvergeom = w3d.RZgeom

w3d.bound0 = periodic
w3d.boundnz = periodic
w3d.boundxy = dirichlet

package("w3d")
generate()

diags = []
fieldDiags = []
solvers = []

if 'electrostatic' in solvertype:
    if w3d.solvergeom == w3d.XYZgeom:
        solver = MultiGrid3D()
    elif w3d.solvergeom == w3d.RZgeom:
        solver = MultiGrid2D()
    solvers.append(solver)
    fieldDiags += [FieldDiagnostic.ElectrostaticFields(
        solver=solver,
        top=top,
        w3d=w3d,
        comm_world=comm_world,
        period=fieldperiod,
        write_dir=diagDir + '/fields/electric'
    )]
if 'magnetostatic' in solvertype:
    solver = MagnetostaticMG()
    solvers.append(solver)
    fieldDiags += [FieldDiagnostic.MagnetostaticFields(
        solver=solver,
        top=top,
        w3d=w3d,
        comm_world=comm_world,
        period=fieldperiod,
        write_dir=diagDir + '/fields/magnetic'
    )]

for solver in solvers:
    registersolver(solver)
    if hasattr(solver, 'mgtol'):
        if type(solver.mgtol) is float:
            solver.mgtol = solvertolerance
        else:
            solver.mgtol = np.array([solvertolerance] * len(solver.mgtol))

if outputFields is True:
    diags += fieldDiags


def generateDist(npart=ptcl_per_step, zemit=dz/5, dv_over_v=0):
    ptclTrans = createKV(
        npart=npart,
        a=0.010,
        b=0.010,
        emitx=4. * 1.e-6,
        emity=4. * 1.e-6
    )

    if dv_over_v > 0:
        vzoffset = np.random.normal(scale=dv_over_v, size=npart)
    else:
        vzoffset = np.reshape(np.zeros(npart), (npart, 1))

    vz = beam_beta * clight * (np.ones(npart) + vzoffset)
    zemit = zemit * np.ones(npart)
    return np.column_stack((ptclTrans, zemit, vz))


def injectelectrons(npart=ptcl_per_step, zoffset=w3d.dz/5):  # emitting surface 1/5th of a cell forward by default
    ptclArray = generateDist(npart=npart, zemit=zoffset, dv_over_v=0.001)
    beam.addparticles(
        x=ptclArray[:, 0],
        y=ptclArray[:, 2],
        z=ptclArray[:, 4],
        vx=ptclArray[:, 1] * ptclArray[:, 5],  # Change to x and y angles to velocities
        vy=ptclArray[:, 3] * ptclArray[:, 5],
        vz=ptclArray[:, 5]
    )


def preloadelectrons():
    for i in range(0, w3d.nz):
        injectelectrons(npart=int(34500/w3d.nz), zoffset=i*dz)

# Install injector
installuserinjection(injectelectrons)

####################################
# Ionization of background gas     #
####################################

ioniz = Ionization(
    stride=100,
    xmin=w3d.xmmin,
    xmax=w3d.xmmax,
    ymin=w3d.ymmin,
    ymax=w3d.ymmax,
    zmin=(w3d.zmmin + w3d.zmmax)/2. - w3d.dz*3,
    zmax=(w3d.zmmin + w3d.zmmax)/2. + w3d.dz*3,
    nx=w3d.nx,
    ny=w3d.ny,
    nz=w3d.nz,
    l_verbose=True
)

target_pressure = 1  # in Pa
target_temp = 273  # in K
target_density = target_pressure / boltzmann / target_temp  # in 1/m^3

if simulateIonization is True:
    # e + H2 -> 2e + H2+
    ioniz.add(
        incident_species=beam,
        emitted_species=[h2plus, emittedelec],
        emitted_energy0=[lambda vi, nnew: h2crosssections.ejectedEnergy(vi, nnew)/10, h2crosssections.ejectedEnergy],  # Array of emission energies (in eV) corresponding to emitted_species
        l_remove_target=False,  # Flag for removing target particle
        cross_section=h2_ioniz_crosssection,
        ndens=target_density
    )

derivqty()  # Sets addition derived parameters (such as beam.vbeam)

##########################
# Injection Controls     #
##########################

# --- Specify injection of the particles
top.inject = 6                       # 2 means space-charge limited injection, 6 is user specified
top.npinject = ptcl_per_step           # Approximate number of particles injected each step
top.ainject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything
top.binject = 0.0008                      # Must be set even for user defined injection, doesn't seem to do anything

############################
# Particle Diagnostics     #
############################

diagP = ParticleDiagnostic(
    period=1,
    top=top,
    w3d=w3d,
    species={species.name: species for species in listofallspecies},
    comm_world=comm_world,
    lparallel_output=False,
    write_dir=diagDir+'/xySlice/'
)
diags += [diagP]


def writeDiagnostics():
    for d in diags:
        d.write()

installafterstep(writeDiagnostics)

package("w3d")
generate()

preloadelectrons()  # load full beam in before first step
step(1)

writeDiagnostics()
diagP.period = particleperiod

stept(args.time)

# Print a bell 3 times to indicate end of run
for i in range(0, 3):
    print(u'\a')
    time.sleep(1)
