#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 09:38:39 2019

@author: coleman
"""

import numpy as np
from scipy import constants,integrate,stats
import matplotlib.pyplot as plt

#Computing integrals of the transverse & longitudinal forces from the 
# D&S model.

#I expect i'll need to refactor this into an OO form for exporting to
# JSPEC, but for now I just want to do the integral

def main():
    
    #Define fundamental constants
    Z       = 40       #Atomic number of the ion in the beam
    r_e     = constants.physical_constants["classical electron radius"]
    m_e     = constants.m_e
    eq      = constants.e
    n_e     = 2e15     # /m^3, electron density
    omega_p = np.sqrt(4*constants.pi*n_e*eq**2/m_e)

    #The scale factor common to all of the force calculations    
    F_const = (-2 * Z**2 * np.pi * n_e)/m_e
    
    #Constants of the beamline
    B = 5           #Tesla, magnetic field strength (1 or 5)
    l_cool  = 10    #Length of the cooling section (Default in Sirepo JSPEC)


    Delta_e = 1E5   # m/s, RMS electron velocity spread, per Ilya's talk
    Delta_trans = 400   # Transverse spread in e velocity (what's a reasonable value?)
    Delta_long = 1e5    # Longitudinal spread in e velocity.
    
    #Directional components and magnitude of ion beam
    V_trans = 4.2e5  # m/s, per Ilya's talk
    V_long  = 1e5    # m/s, per Ilya's talk
    V_ion   = np.sqrt(V_trans**2 + V_long**2)

    #Impact parameters
    p_L   = constants.c * m_e * Delta_trans/(eq*B)
    p_sh  = np.sqrt(V_ion**2 + Delta_e**2)/omega_p
    p_max = min(max(p_sh,(3*Z/n_e)**(1/3)),V_ion*tau)

    #Coulomb log (depends on v_e)
    L_M = np.log(p_max / p_L)
    

    
    #The integral performed is different depending on the impact parameter
    # and the relative velocities of ion and electron
    def integrand_transverse(v_e):
        #This function relies on Z, n_e, and the various velocities
        # to be defined above, and in scope
        
        U = np.sqrt(V_trans**2 + (V_long-v_e)**2)
        
        integrand = (V_trans * (V_trans**2 - 2*(V_long - v_e)**2))/(U**5)
        integrand *= stats.norm.pdf(v_e,0,Delta_e)
        
        return integrand
        
    def integrand_longitudinal(v_e):

        #These are dependent on v_e in the integrand
        gamma   = 107   #Per Ilya's talk for 55MeV electrons
        beta    = np.sqrt(1.0-(1.0/gamma**2))
        #ion time of flight in the cooling section
        tau     = l_cool/(beta*gamma*constants.c)
        v_mean  = beta * constants.c  #Mean of electron gaussian speed distribution

        #Impact parameters
        p_L   = constants.c * m_e * Delta_trans/(eq*B)
        p_sh  = np.sqrt(V_ion**2 + Delta_e**2)/omega_p
        p_max = min(max(p_sh,(3*Z/n_e)**(1/3)),V_ion*tau)

        #Coulomb log (depends on v_e)
        L_M = np.log(p_max / p_L)

        U = np.sqrt(V_trans**2 + (V_long-v_e)**2)
        
        if V_trans <1.0: #Actually at V_trans==0, but lets head off the divergence
            #This actually depends on the longitudinal spread in electron velocity,
            # but we're assuming a large magnetic field so I'm assuming the
            # transverse spread is 0 so Delta_e = Delta_E_longitudinal
            integrand = -V_long * 2 * F_const / (np.sqrt(2*constants.pi)*Delta_e**3)
            #Set this up to integrate over a delta function
            if v_e == V_long:
                integrand *= np.exp(- V_long**2 / 2(Delta_e**2))
            else:
                integrand = 0.0                                                    
                           
        elif V_ion > (100*Delta_long):                      
        
            #When V_ion >> Delta_trans
            integrand = (3 * L_M * V_trans**3 * (V_long-v_e))/(U**5) 
            integrand += 2 * (V_long-v_e) / (U**3)

        else:
            #When V_ion />> Delta_trans
            integrand = (3 * V_trans*(V_long - v_e))/(U**5)        
            integrand *= stats.norm.pdf(v_e,0,Delta_e)
               
        
        
        return integrand
    
        
    #Perform the test integrals
    F_trans = np.float64(F_const) * L_M 
    F_trans *= integrate.quad(lambda x:integrand_transverse(x),-np.inf,np.inf)[0]

    F_long = np.float64(F_const) * L_M 
    F_long *= integrate.quad(lambda x: integrand_longitudinal(x),-np.inf,np.inf)[0]

    print(F_trans)
    print(F_long)
    
     
if __name__ == "__main__":
    main()