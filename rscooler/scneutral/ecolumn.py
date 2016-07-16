"""
From the minimum change from the PXIE LEBT input
IOTA e-column simulation
"""

from warp import *
from egun_like import *
from timedependentvoltage import TimeVoltage
from ionization import *
from extpart import ZCrossingParticles
import sys
#from Secondaries import *

# which geometry to use 2d or 3d
#w3d.solvergeom = w3d.RZgeom
w3d.solvergeom = w3d.XYZgeom


#define some strings that go into the output file
top.pline1     = "e-column"
top.pline2     = "e-column - simplified"
top.runmaker   = "Chong Shik Park (cspark@fnal.gov)"

# --- Invoke setup routine for the plotting
#setup() 
setup(makepsfile=1, cgmlog=1)

# --- Set basic beam parameters
#ibeaminit     =  5e-3   # in A, 5mA
ibeaminit      =  8e-3   # in A, 8mA from HINS
ekininit       =  2.5e6  # in eV, 2.5 MeV HINS proton beam

# Scan parameters
Vbias = -2000.   # Should be negative 
Bsol  = 0.5 

# --- Set input parameters describing the 3d simulation
#top.dt = 1.e-11  # time step
#do faster time steps for testing
#top.dt = 2.5e-10  # time step
#do a much faster time steps for testing
#top.dt = 2e-9 # time step

top.dt = 1.5e-11 # considering the Larmor frequency of electrons in 0.5 T 0.25*(2Pi /wc)

sw = 149794 # for 5 particle per injection 
#sw is calculated from (ibeaminint top.dt / q /  top.npinject)

# relativistic correction
beta   = 0.0728557
fselfb_option = clight*beta

# --- define Species
ions       = Species(type=Hydrogen,charge_state=+1,name='p',     weight = sw)
electrons  = Species(type=Electron, name='e-',                   weight = 2*sw)
h2plus     = Species(type=Dihydrogen, charge_state=+1,name='H2+',weight = 2*sw)
#hplus      = Species(type=Hydrogen, charge_state=+1,name='H+')
#hneutral   = Species(type=Hydrogen, charge_state=0, name='H')
#h2neutral  = Species(type=Dihydrogen, charge_state=0,name='H2')

top.finject[0, ions.jslist[0]]      = 1
top.finject[0, electrons.jslist[0]] = 0
top.finject[0, h2plus.jslist[0]]    = 0
#top.finject[0, hplus.jslist[0]]     = 0
#top.finject[0, hneutral.jslist[0]]  = 0
#top.finject[0, h2neutral.jslist[0]] = 0

# --- starting conditions for the ion and electron beam
top.a0       =    0.01  
top.b0       =    0.01  
top.ap0      =    0.0 
top.bp0      =    0.0
top.vbeam    =    .0e0  # In derivqty(), the beam velocity will be calculated by ekininit anyhow.
top.emit     =    2.e-6 # To define distrbtn, it is likely that we need emittance 
top.ibeam    = ibeaminit
top.ekin     = ekininit   # in eV
top.vthz     =    sqrt(2.*abs(ions.charge)/ions.mass * 250) #convert from eV to m/s
top.lrelativ =    true
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

le = 1.0 # length of the ecolumn length

# --- Set field grid size
w3d.xmmin = -top.prwall
w3d.xmmax = +top.prwall
w3d.ymmin = -top.prwall
w3d.ymmax = +top.prwall
w3d.zmmin =  0.0-0.1
w3d.zmmax  = le +0.1  # 1-m-long ecolumn length

step_ini=int((w3d.zmmax/vz)/top.dt)

if w3d.l4symtry: w3d.xmmin = 0.
if w3d.l2symtry or w3d.l4symtry: w3d.ymmin = 0.

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

# --- Specify injection of the particles
top.npmax     = 500    # To laod initial plasma, this line should be commented out ! (by Dave Grote). Also the number is not corret either.
top.inject    = 1      # 0: no injection, 1: constant current, 2: space-charge limited injection
#top.rinject  = 9999.  # Source radius of curvature # not needed

top.npinject  = 5      # Approximate number of particles injected each step, I reduced it from 1552 to 15

top.linj_efromgrid = true  # Turn on transverse E-fields near emitting surface
top.zinject = w3d.zmmin    # initial z of particle injection?
top.ibpush   = 1           # Specifies type of B-advance; 0 - none, 1 - fast
#top.zinject = 9.51*mm
#top.thetainject = 3.0*deg # tilt the source


w3d.distrbtn = "TE" #Pseudo Thermal Equilibrium
w3d.cylinder = true
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

"""
#---- defin beam with real particle data

w3d.l_inj_user_particles_v = true

xdata,ydata,vxdata,vydata,vzdata=getdatafromtextfile('30.5.lebt1.warp',dims=[5,None])

def createrealbeam():
  if w3d.inj_js == ions.jslist[0]:
#    ii=(random.random(len(xdata))<0.9)
#    w3d.npgrp = sum(ii)
    w3d.npgrp=1552
    offset = [0.00 for i in range(w3d.npgrp)]
    gchange('Setpwork3d')
    w3d.xt[:] = xdata + offset
    w3d.yt[:] = ydata
    w3d.uxt[:] = vxdata
    w3d.uyt[:] = vydata
    w3d.uzt[:] = vzdata

installuserparticlesinjection(createrealbeam)
"""

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

# --- Create the field solver.--> This is the one suggested by Dave
#solver = MultiGrid3D()
#registersolver(solver)

# --- Set up fieldsolver - 7 means the multigrid solver
top.fstype     = 7
f3d.mgtol      = 1.e-1 # Poisson solver tolerance, in volts
f3d.mgparam    =  1.5
f3d.downpasses =  2
f3d.uppasses   =  2

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
package("w3d")
generate()

#################################################################################################################

beampipe        = ZCylinderOut(radius=top.prwall,     zlower=w3d.zmmin, zupper=w3d.zmmax,voltage= 0.0,  xcent=0,ycent=0,zcent=0)
electrode_left  = ZCylinderOut(radius=top.prwall*0.9, zlower=0.0-0.05,  zupper=0.0+0.05, voltage= Vbias,xcent=0,ycent=0)
electrode_right = ZCylinderOut(radius=top.prwall*0.9, zlower=le -0.05,  zupper=le +0.05, voltage= Vbias,xcent=0,ycent=0)

installconductors(beampipe)
installconductors(electrode_left + electrode_right)
addsolenoid(zi=0., zf = le, ri=top.prwall, maxbz= Bsol) 

# --- Recalculate the fields
fieldsolve(-1)

target_density = 3.54e19 #1.e-3 torr

# There was an error in the NERSC in this step.
"""
target_density = zeros((w3d.nx+1, w3d.ny+1, w3d.nz+1))
target_density[:,:,:int(w3d.nz/2)]  = 3.54e19 #1.e-3 torr
target_density[:,:,int(w3d.nz/2):]  = 3.54e19 #1.e-3 torr
"""

# --- setup the charge exchange
ioniz = Ionization(stride=5, zmin=0., zmax=le)
# ionization only in the middle
#####################################################################################################################
#ioniz.add(ions, ions, cross_section=1.0e-21, ndens=1.e22, emitted_energy0=0., emitted_energy_sigma=0., emitted_tag=2)
#####################################################################################################################

#ioniz.add_chargeexchange(ions, ions,
#		cross_section=1.e-20,
#		ndens=1.e17,
#		emitted_energy0=0.,
#		emitted_energy_sigma=0.,
#		emitted_tag=1)

# ------------p ionziation: p + H2 -> p + e + H2+
ioniz.add(ions, emitted_species = [h2plus, electrons], emitted_energy0=0.026, emitted_energy_sigma=0, cross_section=1.5e-21, ndens=target_density)
#ioniz.add(ions, emitted_species = [electrons], emitted_energy0=0.026, emitted_energy_sigma=0, cross_section=1.5e-21, ndens=target_density)

#ioniz.add(ions, emitted_species = [h2plus],    emitted_energy0=0.026, emitted_energy_sigma=0, cross_section=1.5e-21, ndens=target_density)
#ioniz.add(ions, emitted_species = [electrons], emitted_energy0=1,     emitted_energy_sigma=0, cross_section=1.5e-21, ndens=target_density)

# ---------- H- detarchment:  H- + H2 -> H + e + H2
#ioniz.add_detachment(incident_species=ions, target_species=h2neutral, emitted_species=[hneutral, h2neutral, electrons], ndens=target_density)

# ------ electron ionization e + H2 -> H2+ + e + e
#ioniz.add_ionization(incident_species=electrons, target_species=h2neutral, emitted_species=[h2plus, electrons], ndens=1.e17)
ioniz.add_ionization(incident_species=electrons, emitted_species=[h2plus, electrons], emitted_energy0=0.026, emitted_energy_sigma=0, cross_section=1.3e-20, ndens = target_density)
#ioniz.add_ionization(incident_species=electrons, emitted_species=[electrons], emitted_energy0=0.026, emitted_energy_sigma=0, cross_section=1.3e-20, ndens = target_density)

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

save_repetition = 500
plot_repetition = 1000
final_iter = 57483 * 2
iter = 0
dorun = true
while(dorun):
    if(iter % plot_repetition == 0):

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

    if (iter % save_repetition == 0):
        nions = ions.getn()
        xions = ions.getx()
        yions = ions.gety()     
        zions = ions.getz()
        pidions = ions.getpid()

        ions_track = open("protons_%06d.txt" % iter,"w")
        for i in range(nions):
            ions_track.write("%g %15.12f %15.12f %15.12f %g\n" % (i, xions[i], yions[i], zions[i], pidions[i]))
        ions_track.close()

        nelectrons = electrons.getn()
        xelectrons = electrons.getx()
        yelectrons = electrons.gety()
        zelectrons = electrons.getz()
        pidelectrons = electrons.getpid()

        electrons_track = open("electrons_%06d.txt" % iter,"w")
        for i in range(nelectrons):
            electrons_track.write("%g %15.12f %15.12f %15.12f %g\n" % (i, xelectrons[i], yelectrons[i], zelectrons[i], pidelectrons[i]))
        electrons_track.close()

        nh2plus = h2plus.getn()
        xh2plus = h2plus.getx()
        yh2plus = h2plus.gety()
        zh2plus = h2plus.getz()
        pidh2plus = h2plus.getpid()

        h2plus_track = open("h2plus_%06d.txt" % iter,"w")
        for i in range(nh2plus):
            h2plus_track.write("%g %15.12f %15.12f %15.12f %g\n" % (i, xh2plus[i], yh2plus[i], zh2plus[i], pidh2plus[i]))
        h2plus_track.close()

        fdenx = open('denx_%06d.txt' % iter,'w')
        fex   = open('ex_%06d.txt' % iter, 'w')
        for i in range(0, w3d.nx+1):
            iave = (iden[i,iycenter,izcenter]+iden[i,iycenter,izcenter-1]+iden[i,iycenter,izcenter+1])/3.0
            eave = (eden[i,iycenter,izcenter]+eden[i,iycenter,izcenter-1]+eden[i,iycenter,izcenter+1])/3.0
            have = (hden[i,iycenter,izcenter]+hden[i,iycenter,izcenter-1]+hden[i,iycenter,izcenter+1])/3.0
            fdenx.write('%e %e %e %e \n' % (xcoord[i], iave, eave, have))
            fex.write('%e %e\n' % (xcoord[i], ex[i]))
        fdenx.close()
        fex.close()

        fdenz = open('denz_%06d.txt' % iter,'w')
        fez   = open('ez_%06d.txt' % iter, 'w')
        for k in range(1, w3d.nz+1-1):
            iave = (iden[ixcenter,ixcenter, k]+iden[ixcenter,iycenter, k-1]+iden[ixcenter,iycenter, k+1])/3.0
            eave = (eden[ixcenter,ixcenter, k]+eden[ixcenter,iycenter, k-1]+eden[ixcenter,iycenter, k+1])/3.0
            have = (hden[ixcenter,ixcenter, k]+hden[ixcenter,iycenter, k-1]+hden[ixcenter,iycenter, k+1])/3.0
            fdenz.write('%e %e %e %e \n' % (zcoord[k], iave, eave, have))
            fez.write('%e %e\n' % (zcoord[k], ez[k]))
        fdenz.close()
        fez.close()
 
    """
    # I Think this routine causes serious delay 
    edenave = 0.0
    eden  = electrons.get_density()
    for k in range(0, w3d.nz+1):
        edenave = edenave + eden[ixcenter,iycenter,k]
    edentime[iter] = edenave / (w3d.nz+1)
    """

    step(1)
    iter = iter + 1
    if (iter >= final_iter):
        dorun = false
       
# Plot for final state: routine to save ex and ez into files 		
fdenx = open('denx.txt','w') # If I use "a" option, NERSC gives repetition of data
for i in range(0, w3d.nx+1):
    iave = (iden[i,iycenter,izcenter]+iden[i,iycenter,izcenter-1]+iden[i,iycenter,izcenter+1])/3.0
    eave = (eden[i,iycenter,izcenter]+eden[i,iycenter,izcenter-1]+eden[i,iycenter,izcenter+1])/3.0
    have = (hden[i,iycenter,izcenter]+hden[i,iycenter,izcenter-1]+hden[i,iycenter,izcenter+1])/3.0
    fdenx.write('%e %e %e %e \n' % (xcoord[i], iave, eave, have))
fdenx.close()
	
fdenz = open('denz.txt','w')
for k in range(1, w3d.nz+1-1):
    iave = (iden[ixcenter,ixcenter, k]+iden[ixcenter,iycenter, k-1]+iden[ixcenter,iycenter, k+1])/3.0
    eave = (eden[ixcenter,ixcenter, k]+eden[ixcenter,iycenter, k-1]+eden[ixcenter,iycenter, k+1])/3.0
    have = (hden[ixcenter,ixcenter, k]+hden[ixcenter,iycenter, k-1]+hden[ixcenter,iycenter, k+1])/3.0
    fdenz.write('%e %e %e %e \n' % (zcoord[k], iave, eave, have))
fdenz.close()


#time history of electron desntiy 
ptitles(titlet = "Electron density (#/m3) at center", titleb = "Time (s)", titlel =" ")
pla(eden_time, col_time,  color= red)
fma()

fdentime = open('dentime.txt','w')
for i in range(0, len(col_time)):
    fdentime.write('%e %e \n' % (col_time[i],eden_time[i]))
fdentime.close()

# phase-space plot to see any instability
ions.ppzvtheta (view = 3,color=red)
ions.ppzvr(view = 4,color=red)
ions.ppxvx(view = 5,color=red)
ions.ppzvz(view = 6,color=red)
fma()
