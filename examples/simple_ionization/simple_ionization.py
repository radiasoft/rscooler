"""
From the minimum change from the PXIE LEBT input
IOTA e-column simulation
"""

import sys
from warp import *
from warp.run_modes.egun_like import *  # Need full path in module import
from warp.utils.timedependentvoltage import TimeVoltage
from warp.particles.ionization import *
from warp.particles.extpart import ZCrossingParticles
from warp.data_dumping.openpmd_diag import particle_diag
from warp.data_dumping.openpmd_diag import field_diag
import numpy


# which geometry to use 2d or 3d
# w3d.solvergeom = w3d.RZgeom
w3d.solvergeom = w3d.XYZgeom


#define some strings that go into the output file
top.pline1     = "simple_ionization"
top.pline2     = "Demonstration of simple ionization, e- + H2 -> e- + e- + H2+"
top.runmaker   = "James Gerity (jgerity@tamu.edu)"

# --- Invoke setup routine for the plotting
#setup()
setup(makepsfile=1, cgmlog=1)

# --- Set basic beam parameters
ibeaminit = 500e-3   # in A, from Fermilab experiment
ekininit = 114e6  # in eV, from Fermilab experiment

# Scan parameters
Vbias = -100.   # Should be negative
Bsol = 0.0

# --- Set input parameters describing the 3d simulation

w_cyclotron = echarge*Bsol/emass # Larmor frequency
top.dt = 1e-9  # fast time step for testing
print("Timestep is %E s" % top.dt)

top.npinject  = 5.  # Approximate number of particles injected each step
sw = ibeaminit * top.dt / echarge / top.npinject

# relativistic correction
gamma = ekininit / 511.e6
beta = 1. - 1./gamma**2
print("Relativistic factors are gamma=%.5f, beta=%.5f" % (gamma,beta))

# --- define Species
ions = Species(type=Hydrogen, charge_state=+1, name='p', weight=sw)
electrons = Species(type=Electron, name='e-', weight= sw)
h2plus = Species(type=Dihydrogen, charge_state=+1, name='H2+', weight= sw)

# jslist is list of species, top.finject[n,i] is fraction of species i to inject at source n
# Inject only electrons
top.finject[0, ions.jslist[0]]      = 0
top.finject[0, electrons.jslist[0]] = 1
top.finject[0, h2plus.jslist[0]]    = 0

# --- starting conditions for the ion and electron beam
top.a0       = 0.01
top.b0       = 0.01
top.ap0      = 0.0
top.bp0      = 0.0
top.vbeam    = .0e0  # In derivqty(), the beam velocity will be calculated by ekininit anyhow.
top.emit     = 2.e-6 # To define distrbtn, it is likely that we need emittance
top.ibeam    = ibeaminit
top.ekin     = ekininit   # in eV
top.vthz     = 0.0 #convert from eV to m/s
top.lrelativ = true
derivqty()

w3d.l4symtry = false  # four-fold symmetry
w3d.l2symtry = false  # two-fold symmetry

vz=sqrt(2.*abs(ions.charge)*ions.ekin/ions.mass) # ion velocity

# --- Set boundary conditions

# ---   for field solve
w3d.bound0  = neumann #dirichlet
w3d.boundnz = neumann
w3d.boundxy = dirichlet #neumann

# ---   for particles
top.pbound0  = absorb
top.pboundnz = absorb
top.prwall   = 25.e-3 # 25 mm

le = 1.0 # length of the domain

# --- Set field grid size
w3d.xmmin = -top.prwall
w3d.xmmax = +top.prwall
w3d.ymmin = -top.prwall
w3d.ymmax = +top.prwall
w3d.zmmin =  0.0 - 0.1
w3d.zmmax  = le + 0.1

if w3d.l4symtry:
    w3d.xmmin = 0.
if w3d.l2symtry or w3d.l4symtry:
    w3d.ymmin = 0.

# set grid spacing
w3d.dx = (2.*top.prwall)/100.
w3d.dy = (2.*top.prwall)/100.
w3d.dz = 1./100.

# --- Field grid dimensions - nx and ny do not need to be even
w3d.nx    = int((w3d.xmmax - w3d.xmmin)/w3d.dx)
w3d.xmmax = w3d.xmmin + w3d.nx*w3d.dx
w3d.ny    = int((w3d.ymmax - w3d.ymmin)/w3d.dy)
w3d.ymmax = w3d.ymmin + w3d.ny*w3d.dy
w3d.nz    = int((w3d.zmmax - w3d.zmmin)/w3d.dz)
w3d.zmmax = w3d.zmmin + w3d.nz*w3d.dz

print("Grid size is %i x %i x %i" % (w3d.nx,w3d.ny,w3d.nz))

# --- Specify injection of the particles
# --- npmax is the number of simulation particles to create.
top.npmax     = 50000
top.inject    = 1      # 0: no injection, 1: constant current, 2: space-charge limited injection (Child-Langmuir)

top.linj_efromgrid = true  # Turn on transverse E-fields near emitting surface
top.zinject = w3d.zmmin    # initial z of particle injection
top.ibpush   = 1           # Specifies type of B-advance; 0 - none, 1 - fast

w3d.distrbtn = "TE" #Pseudo Thermal Equilibrium
w3d.cylinder = true
# from w3d.F: --- use random numbers to load particles in polar coordinates
w3d.ldprfile = "polar"

top.inject = 1
source_radius = 5.5e-3 # 5.5 mm from MAD

top.ainject = source_radius
top.binject = source_radius
w3d.l_inj_user_particles_v = true

def nonlinearsource():
    if w3d.inj_js == ions.jslist[0]:
        np = top.npinject
        r = source_radius*random.random(np)
        theta = 2.*pi*random.random(np)
        x = r*cos(theta)
        y = r*sin(theta)
        w3d.npgrp = np
        gchange('Setpwork3d')
        w3d.xt[:] = x
        w3d.yt[:] = y
        w3d.uxt[:] = 0.
        w3d.uyt[:] = 0.
        w3d.uzt[:] = vz

installuserparticlesinjection(nonlinearsource)

# --- Select plot intervals, etc.
top.nhist = 1 # Save history data every time step
top.itplfreq[0:4]=[0,1000000,25,0] # Make plots every 25 time steps
top.itmomnts[0:4]=[0,1000000,top.nhist,0] # Calculate moments every step

# --- Save time histories of various quantities versus z.
top.lhcurrz  = true
top.lhrrmsz  = true
top.lhxrmsz  = true
top.lhyrmsz  = true
top.lhepsnxz = true
top.lhepsnyz = true
top.lhvzrmsz = true

# --- Set up fieldsolver - 7 means the multigrid solver
top.fstype     = 7
f3d.mgtol      = 1.e-1 # Poisson solver tolerance, in volts
f3d.mgparam    =  1.5
f3d.downpasses =  2
f3d.uppasses   =  2

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
package("w3d")
generate()
h5output = 'diags/hdf5/'

def cleanupPrevious(outputDirectory = 'diags/hdf5/'):
    if os.path.exists(outputDirectory):
        files = os.listdir(outputDirectory)
        for file in files:
            if file.endswith('.h5'):
                os.remove(os.path.join(outputDirectory,file))

cleanupPrevious()

#Version of Warp currently installed on Container doesn't handle species weights correctly - 4/21/16
diagP = particle_diag.ParticleDiagnostic( period=5, top=top, w3d=w3d,
        species= { species.name : species for species in listofallspecies },
        particle_data=["position","momentum"],
        comm_world=comm_world, lparallel_output=True, )

installafterstep( diagP.write )

#################################################################################################################

beampipe        = ZCylinderOut(radius=top.prwall,     zlower=w3d.zmmin, zupper=w3d.zmmax,voltage= 0.0,  xcent=0,ycent=0,zcent=0)
electrode_left  = ZCylinderOut(radius=top.prwall*0.9, zlower=0.0-0.05,  zupper=0.0+0.05, voltage= Vbias,xcent=0,ycent=0)
electrode_right = ZCylinderOut(radius=top.prwall*0.9, zlower=le -0.05,  zupper=le +0.05, voltage= Vbias,xcent=0,ycent=0)

installconductors(beampipe)
installconductors(electrode_left + electrode_right)

# --- Recalculate the fields
fieldsolve(-1)

target_density = 3.54e19 #1.e-3 torr

# --- setup the charge exchange
# From ionization.py:
# stride=100: The stride in the particle loops, controlling how many particles
#               can under go a collision each time step.

ioniz = Ionization(stride=top.npmax / 1000, zmin=0., zmax=le)

# ------------ e + H2 -> e + e + H2+
ioniz.add(incident_species=electrons,
          emitted_species=[h2plus, electrons],
          emitted_energy0=[0.026, 0.026], # Array of emission energies corresponding to emitted_species
          emitted_energy_sigma=0,
          l_remove_target=False, # Flag for removing target particle
          cross_section=1.5e-21, # NIST total ionization cross-section data suggests this is probably ~ 0.01
          ndens=target_density)

ixcenter = int(w3d.nx/2)
iycenter = int(w3d.ny/2)
izcenter = int(w3d.nz/2)

xcoord = zeros(w3d.nx+1,float)
for k in range(0, w3d.nx+1):
    xcoord[k] = w3d.xmmin + k * w3d.dx

zcoord = zeros(w3d.nz+1,float)
for k in range(0, w3d.nz+1):
    zcoord[k] = w3d.zmmin + k * w3d.dz

eden_time = []
col_time = []

def doRun():
    num_steps = 500
    run_time = 1.e-6
    final_iter = run_time / top.dt
    iter = 0
    while(iter < final_iter):
        ions.ppzx(color=red, titles = 0, view=9, pplimits=(w3d.zmmin, w3d.zmmax, -top.prwall, top.prwall))
        pfzx(fill=1,filled=0, plotsg=0,titles = 0, cond=0,view=9)

        h2plus.ppzx(color=blue, titles = 0,msize=100, view=10)
        electrons.ppzx(color=green, titles = 0, msize=100, view=10, pplimits=(w3d.zmmin, w3d.zmmax, -top.prwall, top.prwall))
        fma()

        iden  = ions.get_density()
        eden  = electrons.get_density()
        hden  = h2plus.get_density()

        ptitles(titlet = "Density (#/m3) at x=y=0", titleb = "Z (m)", titlel =" ")
        limits(w3d.zmmin, w3d.zmmax)
        pla(iden[ixcenter,iycenter,0:], zcoord, color = red)
        pla(eden[ixcenter,iycenter,0:], zcoord, color = green)
        pla(hden[ixcenter,iycenter,0:], zcoord, color = blue)
        fma()

        ptitles(titlet = "Density (#/m3) at ecolumn center", titleb = "X (m)", titlel =" ")
        limits(w3d.xmmin, w3d.xmmax)
        pla(iden[0:,iycenter,izcenter], xcoord, color = red)
        pla(eden[0:,iycenter,izcenter], xcoord, color = green)
        pla(hden[0:,iycenter,izcenter], xcoord, color = blue)
        fma()

        # longitudinal electric fields
        ez = getselfe(comp="z", ix = ixcenter,  iy = iycenter)
        ptitles(titlet = "Self electric fields (V/m) along the beam", titleb = "Z (m)",titlel =" ")
        limits(w3d.zmmin, w3d.zmmax)
        pla(ez, zcoord,  color= red)
        fma()

        # transverse electric fields
        ex = getselfe(comp="x", iy = iycenter,  iz = izcenter)
        ptitles(titlet = "Self electric fields (V/m) along the beam", titleb = "X (m)",titlel =" ")
        limits(w3d.xmmin, w3d.xmmax)
        pla(ex, xcoord,  color= red)
        fma()

        eden_time.append( (eden[ixcenter,iycenter,izcenter] + eden[ixcenter,iycenter,izcenter-1] + eden[ixcenter,iycenter,izcenter+1])/3. )
        col_time.append( top.time )

        step(num_steps)
        iter = iter + num_steps

    #time history of electron desntiy
    ptitles(titlet = "Electron density (#/m3) at center", titleb = "Time (s)", titlel =" ")
    pla(eden_time, col_time,  color= red)
    fma()
