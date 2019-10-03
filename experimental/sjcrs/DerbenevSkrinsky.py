#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 09:38:39 2019

@author: coleman
"""

import numpy as np
from scipy import constants,integrate,stats
import matplotlib.pyplot as plt
import seaborn as sns

#Computing integrals of the transverse & longitudinal forces from the 
# D&S model.

class ebunch:
    def __init__(self):
        #Define fundamental constants
        self.m_e     = constants.m_e
        self.eq      = constants.e
        self.n_e     = 2e15     # /m^3, electron density
        self.omega_p = np.sqrt(4*constants.pi*self.n_e*self.eq**2/self.m_e)

        self.v_e = 0.0 #Provide a space for v_e to mutate

        #The scale factor common to all of the force calculations    
        self.F_const = np.float64((-2 * np.pi * self.n_e)/self.m_e)
    
        #Constants of the beamline
        self.B = 5           #Tesla, magnetic field strength (1 or 5)
        self.l_cool  = 10    #Length of the cooling section (Default in Sirepo JSPEC)

        #These are RMS spreads of e velocity
        self.Delta_trans = 400   # Transverse spread(what's a reasonable value?)
        self.Delta_long = 1e5    # Longitudinal spread
        self.Delta_e = np.sqrt(self.Delta_long**2 + self.Delta_trans**2)
    
        #these are v_e dependent, but I'll keep them constant for now
        self.gamma   = 107   #Per Ilya's talk for 55MeV electrons
        self.beta    = np.sqrt(1.0 - (1.0 / self.gamma**2))
        #ion time of flight in the cooling section
        self.tau     = self.l_cool/(self.beta * self.gamma * constants.c)
        self.v_mean  = self.beta * constants.c  #Mean of electron gaussian speed distribution


class ionbunch:
    def __init__(self):
        self.Z       = 79       #Atomic number of the ion in the beam (gold=79)    

        #Directional components and magnitude of ion beam
        self.V_trans = 4.2e5  # m/s, per Ilya's talk
        self.V_long  = 1e5    # m/s, per Ilya's talk
        self.V_ion   = np.sqrt(self.V_trans**2 + self.V_long**2)

    #Impact parameters
    def ImpactParameters(self,EB): #EB is an electronbunch object    

        #The value of the RMS spreads from EB determines p_max
        p_L   = constants.c * constants.m_e * EB.Delta_trans / (EB.eq*EB.B)
        p_sh  = np.sqrt(self.V_ion**2 + EB.Delta_e**2) / EB.omega_p
        p_max = min(max(p_sh,(3*self.Z/EB.n_e)**(1/3)),self.V_ion*EB.tau)

        #Coulomb log (depends on v_e)
        L_M = np.log(p_max / p_L)
        return L_M

    #The integral performed is different depending on the impact parameter
    # and the relative velocities of ion and electron
    def integrand_transverse(self,v_e,EB):
        
        U = np.sqrt(self.V_trans**2 + (self.V_long-v_e)**2)
        
        EB.v_e = v_e
        
        L_M = self.ImpactParameters(EB)
        
        integrand = (self.V_trans * (self.V_trans**2 - 2*(self.V_long - v_e)**2))/(U**5)
        integrand *= stats.norm.pdf(v_e,0,EB.Delta_e)
        
        return L_M * integrand
        
    def integrand_longitudinal(self,v_e,EB):

        U = np.sqrt(self.V_trans**2 + (self.V_long-v_e)**2)
        
        EB.v_e = v_e        
        
        L_M = self.ImpactParameters(EB)
        
        if self.V_trans < 1.0: 
            #Actually at V_trans==0, but lets head off the divergence
            # by cutting of near 0
            
            integrand = -self.V_long * 2 * EB.F_const 
            integrand /= np.sqrt(2 * constants.pi) * EB.Delta_e**3
            #Set this up to integrate over a delta function
            if v_e == self.V_long:
                integrand *= np.exp(- self.V_long**2 /(2 * EB.Delta_e**2))
            else:
                integrand = 0.0                                                    
                           
        elif self.V_ion > (100.0 * EB.Delta_long):                      
            #When V_ion >> Delta_long
            integrand = (3.0 * self.V_trans**2 * (self.V_long-v_e))
            integrand /= (U**5) 
            integrand += 2.0 * (self.V_long-v_e) / (U**3)
            integrand *= stats.norm.pdf(v_e, 0, EB.Delta_e)

        else:
            #When V_ion />> Delta_long
            integrand = (3.0 * self.V_trans**2 * (self.V_long - v_e))
            integrand /= (U**5)        
            integrand *= stats.norm.pdf(v_e, 0, EB.Delta_e)            
        
        return L_M * integrand

    def Force_transverse(self,EB):
                
        F_trans = EB.F_const * self.Z**3
        #integration only happens over the first variable
        # tuple returned is (value,error). Just grab the value
        F_trans *= integrate.quad(lambda x:self.integrand_transverse(x,EB),-np.inf,np.inf)[0]

        return F_trans        

    def Force_longitudinal(self,EB):
        
        F_long = EB.F_const * self.Z**3
        F_long *= integrate.quad(lambda x: self.integrand_longitudinal(x,EB),-np.inf,np.inf)[0]
        
        return F_long

        
def main():

    EB = ebunch()
    #Change any initial conditions if you want
    
    vec_trans = []
    vec_long = []
    
    for i,x in enumerate(np.linspace(0,6e5,20)):
        
        IB = ionbunch()
        IB.v_long = x
        IB.v_trans = 0.0
        IB.Z = 1
    
        F_trans = IB.Force_transverse(EB)       
        F_long = IB.Force_longitudinal(EB)

        if i < 5:
            print(F_trans)
            print(F_long)
    
        vec_trans.append(F_trans)
        vec_long.append(F_long)
     
    plt.subplot(211)
    plt.plot(np.linspace(0,6e5,20),vec_long)
    plt.xlabel(r"$V_{ion,\parallel}$(m/s)",fontsize=15)
    plt.ylabel(r"$-F_{\parallel}$ (eV/m)",fontsize=15)
    plt.subplot(212)
    plt.plot(np.linspace(0,6e5,20),vec_trans)
    plt.xlabel(r"$V_{ion,\bot}$(m/s)",fontsize=15)
    plt.ylabel(r"$-F_{\bot}$ (eV/m)",fontsize=15)
    plt.tight_layout()
    plt.show()
        
if __name__ == "__main__":
    main()