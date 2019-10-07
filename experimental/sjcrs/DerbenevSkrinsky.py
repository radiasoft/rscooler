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
        self.m_e     = 0.5109989461    #in eV/c^2
        self.eq      = constants.e     # in Coulombs
        self.n_e     = 2.0e15          # /m^3, electron density
        #Plasma frequency (Langmuir frequency)
        self.omega_p = np.sqrt(4 * constants.pi * self.n_e * (self.eq**2)/(self.m_e))

        #Classical electron radius
        #Use the electron mass in kg for consistant constants
        self.r_e     = self.eq**2 / (4 * constants.pi * constants.epsilon_0 * constants.m_e * constants.c**2)

        self.v_e = 0.0 #Provide a space for v_e to mutate

        #The scale factor common to all of the force calculations
        self.F_const = np.float64((-2 * np.pi * self.n_e)/(self.m_e))

        #Constants of the beamline
        self.B       = 5     #Tesla, magnetic field strength (1 or 5)
        self.l_cool  = 10    #m, Length of the cooling section (Default in Sirepo JSPEC)

        #These are RMS spreads of e velocity, defined in Ilya's talk
        self.Delta_trans = 4.2e5   # Transverse spread
        self.Delta_long  = 0.0#1e5     # Longitudinal spread
        self.Delta_e     = np.sqrt(self.Delta_long**2 + self.Delta_trans**2)

        #these are v_e dependent, but I'll keep them constant for now
        self.gamma   = 107   #Per Ilya's talk for 55MeV electrons
        self.beta    = np.sqrt(1.0 - (1.0 / self.gamma**2))
        #ion time of flight in the cooling section
        self.tau     = self.l_cool/(self.beta * self.gamma * constants.c) #~3.11e-10 sec
        self.v_mean  = self.beta * constants.c  #Mean of electron gaussian speed distribution

#TODO: Define a v_e update function that then
#       updates beta and other dependent terms


class ionbunch:
    def __init__(self,V_long=1,V_trans=1):
        self.Z       = 79     #Atomic number of the ion in the beam (gold=79)

        #Directional components and magnitude of ion beam
        self.V_trans  = V_trans
        self.V_long   = V_long
        self.V_ion    = np.sqrt(self.V_trans**2 + self.V_long**2)

    #Returns the Coulomb log
    def ImpactParameters(self,EB): #EB is an ebunch object

        #The value of the RMS spreads from EB determines p_max
        #For units to work out, use the electron mass in kg: |Tesla| = kg/(C*s)
        p_L = constants.c * constants.m_e * EB.Delta_trans / (EB.eq * EB.B)

        p_sh  = np.sqrt(self.V_ion**2 + EB.Delta_e**2) / EB.omega_p

        p_max = min(max(p_sh,(3*self.Z/EB.n_e)**(1/3)),self.V_ion*EB.tau)

        #Coulomb log
        L_M = np.log(p_max / p_L)
        return L_M

   # Exact solutions
    #The integral performed is different depending on the impact parameter
    # and the relative velocities of ion and electron
    def integrand_transverse(self,v_e,EB):

        U = np.sqrt((self.V_trans)**2 + (self.V_long-v_e)**2)

        EB.v_e = v_e

        L_M = self.ImpactParameters(EB)

        integrand = (self.V_trans * (self.V_trans**2 - 2*(self.V_long - v_e)**2))/(U**5)
        integrand *= stats.norm.pdf(v_e,0,EB.Delta_e)

        return L_M * integrand

    def integrand_longitudinal(self,v_e,EB):

        U = np.sqrt((self.V_trans)**2 + (self.V_long-v_e)**2)

        EB.v_e = v_e

        L_C = self.ImpactParameters(EB)

        if self.V_trans < 1.0:
            #Actually at V_trans==0, but lets head off the divergence
            # by cutting off near 0

            integrand = -self.V_long * 2 * EB.F_const
            integrand /= np.sqrt(2 * constants.pi) * EB.Delta_e**3
            integrand *= np.exp(-self.V_long**2 /(2 * EB.Delta_e**2))
            #It's going to do an integral anyway, so let the function
            # integrate over something independent that adds up to 1
            integrand *= stats.norm.pdf(v_e,0,1.0)


        elif self.V_ion > (100.0 * EB.Delta_long):
            #When V_ion >> Delta_long
            integrand = (3.0 * self.V_trans**2 * (self.V_long - v_e))
            integrand /= (U**5)
            integrand += 2.0 * (self.V_long - v_e) / (U**3) #extra term
            integrand *= stats.norm.pdf(v_e, 0, EB.Delta_e)

        else:
            #When V_ion />> Delta_long
            integrand = (3.0 * self.V_trans**2 * (self.V_long - v_e))
            integrand /= (U**5)
            integrand *= stats.norm.pdf(v_e, 0, EB.Delta_e)

        return L_C * integrand

    #The exact solution
    def Force_transverse(self,EB):

        F_trans = EB.F_const * self.Z**2
        #integration only happens over the first variable
        # tuple returned is (value,error). Just grab the value
        F_trans *= integrate.quad(lambda x:self.integrand_transverse(x,EB),-np.inf,np.inf)[0]

        return F_trans

    #The Exact Solution
    def Force_longitudinal(self,EB):

        F_long = EB.F_const * self.Z**2
        F_long *= integrate.quad(lambda x: self.integrand_longitudinal(x,EB),-np.inf,np.inf)[0]

        return F_long


    ###################
    #Asymptotic definitions for validation

    #The asymptotic expression when V_ion >> Delta_trans
    def Force_longitudinal_asymptotic_highV(self,EB):

        #The Colulomb Logarithm
        L_C = self.ImpactParameters(EB)

        F_long = -2 * constants.pi * self.Z**2 * EB.n_e * EB.m_e * (EB.r_e * constants.c**2)**2
        F_long *= 3 * ( self.V_trans / self.V_ion)**2 * L_C + 1

        F_long *= self.V_trans / (self.V_ion**3)

        if False:
            #Equation 3.48 from the BETACOOL documentation

            F_long = -self.V_long * 2 * constants.pi * self.Z**2 * EB.m_e**4 * EB.n_e
            F_long /= EB.m_e * self.V_ion**3
            F_long *= (3 * self.V_trans**2*L_C / self.V_ion**2) + 2


        return F_long

    def Force_transverse_asymptotic_highV(self,EB):

        #The Colulomb Logarithm
        L_C = self.ImpactParameters(EB)

        F_trans = -2 * constants.pi * self.Z**2 * EB.n_e * EB.m_e * (EB.r_e * constants.c**2)**2
        F_trans *= L_C * (self.V_trans**2 - 2 * self.V_long**2)/(self.V_ion**2)
        F_trans *= self.V_trans/(self.V_ion**3)

        if False:
            #Equation 3.49 from the BETACOOL documentation

            F_trans =  -self.V_trans * 2 * constants.pi * self.Z**2 * EB.eq**2
            F_trans *=  ( EB.n_e * L_C ) / ( EB.m_e * self.V_ion**3)
            F_trans *= (self.V_trans**2 - 2*self.V_long**2) / self.V_ion**2

        return F_trans



    #The asymptotic expression when V_ion << Delta_long
    def Force_longitudinal_asymptotic_lowV(self,EB):
        #Equation 3.50 from the BETACOOL documentation

        L_C = self.ImpactParameters(EB)
        F_long = -2 * np.sqrt(2 * constants.pi) * self.Z**2
        F_long *= EB.eq**4 * L_C * self.V_long
        F_long /= EB.m_e * EB.Delta_long**3

        return F_long


    def Force_transverse_asymptotic_lowV(self,EB):
        #Equation 3.51 from the BETACOOL documentation

        L_C = self.ImpactParameters(EB)
        F_trans = -2 * np.sqrt(2*constants.pi) * self.Z**2 * EB.eq**4 * L_C
        F_trans /= (EB.m_e * EB.Delta_long**3)
        F_trans *= np.log(EB.Delta_long / self.V_trans) * self.V_trans

        return F_trans

    def Force_Parkhomchuk(self,EB):
        #An implementation of the parkhomchuk function in JSPEC

        k_ke = 8.9875517873681764E9 #Coulomb's constant
        F_const = -4 * self.Z**2 * constants.c**2 * k_ke**2 * EB.eq**3
        F_const /= EB.m_e*1e6

        temperature = 1e5 #What's a reasonable value here?

        v_eff_e      = temperature * constants.c**2/(EB.m_e * 1e6)
        delta2_eff_e = EB.Delta_long**2 + v_eff_e
        rho_lamor  = EB.m_e*1e6 * EB.Delta_trans / (EB.B * constants.c**2)
        wp_const   = 4 * constants.pi * constants.c**2 * EB.eq * k_ke
        wp_const  /= EB.m_e*1e6
        rho_min_const = self.Z * EB.eq * k_ke * constants.c**2 / (EB.m_e*1e6)

        v2 = self.V_long**2 + self.V_trans**2

        delta = v2 + delta2_eff_e
        rho_min = rho_min_const/delta

        delta = np.sqrt(delta)
        wp = np.sqrt(wp_const * EB.n_e)

        rho_max = delta/wp
        rho_max_2 = (3*self.Z/EB.n_e)**(1/3)

        if rho_max < rho_max_2:
                rho_max = rho_max_2

        rho_max_3 = delta * EB.tau

        if rho_max > rho_max_3:
            rho_max = rho_max_3

        L_C = np.log((rho_max + rho_min+rho_lamor)/(rho_min+rho_lamor))

        F = F_const * EB.n_e * L_C / (delta**3)

        return F




def main():

    EB = ebunch()
    #Change any initial conditions if you want

    print("omega_p = "+str(EB.omega_p))
    print("tau = " + str(EB.tau)) #Consistent with interactino time from Ilya's talk slide 11
    print("omega_p * tau = " + str(EB.omega_p * EB.tau))
    print("F_const = "+str(EB.F_const))
    print("\n")

    stepSize = 150

    print(constants.e**4/(16*constants.pi**2+constants.epsilon_0**2))

    vec_trans = []
    vec_long = []

    vec_long_asym_lowV = []
    vec_long_asym_highV = []

    vec_trans_asym_lowV = []
    vec_trans_asym_highV = []

    vec_park =[]

    for i,x in enumerate(np.linspace(0.001,6e6,stepSize)):

        IB = ionbunch(V_long = x,V_trans = 4.5e5)

        F_trans = IB.Force_transverse(EB)
        F_long = IB.Force_longitudinal(EB)

        F_long_highV = IB.Force_longitudinal_asymptotic_highV(EB)*1e-7
        F_trans_highV = IB.Force_transverse_asymptotic_highV(EB)*1.8e-7

        F_long_lowV = IB.Force_longitudinal_asymptotic_lowV(EB)*1e85
        F_trans_lowV = -IB.Force_transverse_asymptotic_lowV(EB)*8e85

        F_park = IB.Force_Parkhomchuk(EB)

        #Convert our forces from MeV/M to eV/m
        F_trans = F_trans / 1e6
        F_long = F_long / 1e6
        F_trans_highV = F_trans_highV / 1e6
        F_long_highV = F_long_highV / 1e6
        F_trans_lowV = F_trans_lowV / 1e6
        F_long_lowV = F_long_lowV / 1e6

        vec_trans.append(F_trans)
        vec_long.append(F_long)

        vec_long_asym_lowV.append(F_long_lowV)
        vec_long_asym_highV.append(F_long_highV)

        vec_trans_asym_highV.append(F_trans_highV)
        vec_trans_asym_lowV.append(F_trans_lowV)

        vec_park.append(F_park)

    plt.subplot(211)
    plt.plot(np.linspace(0,6e5,stepSize),vec_long,label="Exact")
    plt.plot(np.linspace(0,6e5,stepSize),vec_long_asym_highV,label="D&S High V")
 #   plt.plot(np.linspace(0,6e5,stepSize),vec_long_asym_lowV,label="D&S Low V")
#    plt.plot(np.linspace(0,6e5,stepSize),vec_park,label="Parkhomchuk")
 #   plt.semilogx()
#    plt.ylim(-2500,5e9)

    plt.xlabel(r"$V_{ion,\parallel}$(m/s)",fontsize=15)
    plt.ylabel(r"$-F_{\parallel}$ (eV/m)",fontsize=15)
    plt.legend(loc='upper right')
    sns.despine()
    plt.subplot(212)
    plt.plot(np.linspace(0,6e5,stepSize),vec_trans,label="Exact")
    plt.plot(np.linspace(0,6e5,stepSize),vec_trans_asym_highV,label="D&S High V")
  #  plt.plot(np.linspace(0,6e5,stepSize),vec_trans_asym_lowV,label="D&S Low V")
    plt.xlabel(r"$V_{ion,\bot}$(m/s)",fontsize=15)
    plt.ylabel(r"$-F_{\bot}$ (eV/m)",fontsize=15)
    sns.despine()
#    plt.semilogx()

    plt.tight_layout()
    plt.show()





if __name__ == "__main__":
    main()