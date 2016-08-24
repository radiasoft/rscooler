"""
    A example demonstrating ionization of background H2 gas by an electron beam
"""
from __future__ import division
import numpy as np
import shutil
from shutil import os

from warp import *
from warp.data_dumping.openpmd_diag.particle_diag import ParticleDiagnostic
from warp.data_dumping.openpmd_diag import field_diag

from rswarp.diagnostics import FieldDiagnostic
from rswarp.utilities.beam_distributions import createKV
from rswarp.utilities.ionization import Ionization

from warp.diagnostics import gistdummy as gist

solvertype = []
solvertype += ['EM3D']
# solvertype += ['magnetostatic']
# solvertype += ['electrostatic']
outputFields = True
# outputFields = False
fieldperiod = 10  # number of steps between outputting fields
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

beam_gamma = 116e3/511e3 + 1
beam_beta = np.sqrt(1-1/beam_gamma**2)

top.lrelativ = True
top.relativity = 1

sw = 1000

# beam = Species(type=Electron, name='e-', fselfb=beam_beta * clight, weight=sw)
beam = Species(type=Electron, name='e-', weight=sw)
h2plus = Species(type=Dihydrogen, charge_state=+1, name='H2+', weight=1)
emittedelec = Species(type=Electron, name='emitted e-', weight=sw)

beam.ibeam = 1e-6

################################
# 3D Simulation Parameters     #
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
ptcl_per_step = int(beam.ibeam * top.dt / echarge / sw)  # number of particles to inject on each step

top.ibpush = 2  # 0:off, 1:fast, 2:accurate

# --- Other injection variables
w3d.l_inj_exact = True
w3d.l_inj_area = False

# w3d.solvergeom = w3d.XYZgeom
w3d.solvergeom = w3d.RZgeom

# w3d.bound0 = dirichlet
# w3d.boundnz = dirichlet
# w3d.bound0 = openbc
# w3d.boundnz = openbc
w3d.bound0 = periodic
w3d.boundnz = periodic
# w3d.boundxy = neumann
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
    fieldDiags += [FieldDiagnostic.ElectrostaticFields(solver=solver, top=top, w3d=w3d, period=fieldperiod)]
if 'magnetostatic' in solvertype:
    solver = MagnetostaticMG()
    solvers.append(solver)
    fieldDiags += [FieldDiagnostic.MagnetostaticFields(solver=solver, top=top, w3d=w3d, period=fieldperiod)]
if 'EM3D' in solvertype:
    solver = EM3D()
    solvers.append(solver)
    fieldDiags += [field_diag.FieldDiagnostic(em=solver, top=top, w3d=w3d, period=fieldperiod, write_dir=diagDir)]

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


def injectelectrons(npart=ptcl_per_step, zoffset=dz/5):  # emitting surface 1/5th of a cell forward by default
    ptclArray = generateDist(npart=npart, zemit=zoffset, dv_over_v=0.001)
    # Change to x and y angles to velocities
    beam.addparticles(
        x=ptclArray[:, 0],
        y=ptclArray[:, 2],
        z=ptclArray[:, 4],
        vx=ptclArray[:, 1] * ptclArray[:, 5],
        vy=ptclArray[:, 3] * ptclArray[:, 5],
        vz=ptclArray[:, 5]
    )


def preloadelectrons():
    for i in range(0, w3d.nz):
        injectelectrons(npart=int(34500/w3d.nz),zoffset=i*dz)


# Install injector
installuserinjection(injectelectrons)

####################################
# Ionization of background gas     #
####################################

ioniz = Ionization( stride=100,
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
beam_ke = (beam_gamma-1) * 511e3  # in eV

# ------------ e + H2 -> 2e + H2+

bohr_radius = 5.29177e-11  # Bohr Radius
I = 15.42593  # Threshold ionization energy (in eV), from the NIST Standard Reference Database (via NIST Chemistry WebBook)
R = 13.60569  # Rydberg energy (in eV)
S = 4 * np.pi * bohr_radius**2 * (R/I)**2
fitparametern = 2.4  # species-dependent fitting parameter

def fitF(t):
    # Parameters for H2 specific to this fit
    a1 = 0.74
    a2 = 0.87
    a3 = -0.60
    return 1/t * (a1 * np.log(t) + a2 + 1/t * a3)


def fitf1(w, t, n):
    return 1/(w+1)**n + 1/(t-w)**n - 1/((w+1) * (t-w))**(n/2)


def normalizedKineticEnergy(vi=None):
    """
    Compute the normalized kinetic energy t = T/I given an input velocity
    """
    gamma_in = 1. / np.sqrt(1 - (vi/clight)**2)
    T = (gamma_in - 1) * emass * clight**2 / jperev  # kinetic energy (in eV) of incident electron
    t = T / I  # normalized kinetic energy
    return t

def h2_ioniz_diff_crosssection(vi=None, ke_emitted=None):
    """
    Compute the DIFFERENTIAL cross-section for impact ionization of H2 by e-, per David's formula
    vi - incident electron velocity in m/s; this is passed in from warp as vxi=uxi*gaminvi etc.
    ke_emitted - kinetic energy (in eV) of emitted secondary electron
    """
    t = normalizedKineticEnergy(vi)
    w = ke_emitted / I
    F = fitF
    f1 = fitf1
    sigma = S * F(t) * f1(w, t) / I
    return sigma


def h2_ioniz_crosssection(vi=None):
    """
    Compute the TOTAL cross-section for impact ionization of H2 by e-, per David's formula
    vi - incident electron velocity in m/s; this is passed in from warp as vxi=uxi*gaminvi etc.
    """
    t = normalizedKineticEnergy(vi)
    n = fitparametern
    F = fitF

    def g1(t, n):
        return (1 - t**(1-n)) / (n-1) - (2 / (t+1))**(n/2) * (1 - t**(1 - n/2)) / (n-2)

    sigma = S * F(t) * g1(t, n)
    return sigma

if simulateIonization is True:
    # e + H2 -> 2e + H2+
    ioniz.add(incident_species=beam,
              emitted_species=[h2plus, emittedelec],
              emitted_energy0=[1, 1], # Array of emission energies (in eV) corresponding to emitted_species
              emitted_energy_sigma=[0, 0.1],
              l_remove_target=False, # Flag for removing target particle
              cross_section=h2_ioniz_crosssection,
              ndens=target_density)

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
    write_dir=diagDir+'/xySlice/')
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

stept(0.5e-6)
