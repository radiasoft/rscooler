def main(): 
  
  Tsim = 2.0e-10  # s, in the beam frame 
  Nstep = 1000/2  # number of timesteps 
  
  Nconds_z = 4001  #3001 
  
  Z = 1 
  #Bz = 5.0  # T (working in MKS units)
  n_e = 2.0e+15  # m^-3, in the beam frame 
  Delta_e_par = 1.0e+5  # m/s 
  #Delta_e_tr = 4.2e+5  # m/s 
  
  r1 = np.power(3./(4.*np.pi*n_e), 1./3.)  # expectation value of the distance to the nearest electron, assuming isotropic const density 
  
  frac_delta_e_par = 0.01 #0.75  # ratio of the ion velocity to the thermal electron velocity (abs value) 
  # The fact that Delta_e_par is the rms thermal velocity is not used here; it really is only setting the vel. scale 
  V_i_par = frac_delta_e_par *Delta_e_par  # m/s 
  #vz_ini = -V_i_par  # longitudinally cold electrons, ion stationary at the coord system origin 
  
  q_el = -1.6021766208e-19  # Coulombs 
  m_el = 9.10938356e-31  # kg 
  c0 = 299792458.  # m/s 
  r_cl_el = 2.8179403227e-15  # m, classical radius of electron 
		
  Zrcc = Z *r_cl_el *c0 *c0  #  Z * 253.2638458397924 m^3 / s^2 == Z*e*e /([4*pi*eps0]*m_e) 
		
  #w_L = np.abs(q_el*Bz) /m_el  # magnitude of gyration frequency in MKS 
  #T_L = 2.*np.pi /w_L 
  #r_L = Delta_e_tr /w_L  # a "typical" gyro-radius
  
  Nzr_lo = -2.2  # used for setting the lower boundary of the z domain 
  Nzr_hi =  2.4  # upper boundary 
  z_ini  = np.linspace(Nzr_lo *r1, Nzr_hi *r1, Nconds_z) 
  dz = z_ini[1] -z_ini[0] 
  F_time_aved = 0.0 *z_ini 
  
  Nsamp_D = 50 #30 
  dD = 0.1 *r1 
  
  times_r1 = 0.5 +np.arange(Nsamp_D) 
  D = times_r1 *dD 
  
  F_line_integr = 0.0 *D 
  #T_lin = 2. *np.pi *np.sqrt( D**3 /Zrcc ) 
  
  print ' ' 
  print "Z = ", Z 
  print 'n_e (beam fr.) = ', n_e, 'm^-3' 
  print "r1 == r_nearest (in 3D) = ", r1, "m" 
  
  print 'V_i_par = ', V_i_par, 'm/s, Tsim = ', Tsim, 's'  
  print 'V_i_par *Tsim = ', V_i_par *Tsim, 'm,  V_i_par *Tsim /r1 = ', V_i_par *Tsim /r1
  print 'Nstep = ', Nstep, ', dt = ', Tsim /np.float64(Nstep), 's'
  print 'V_i_par *dt = ', V_i_par *Tsim /np.float64(Nstep), 'm,  dz / (V_i_par *dt)  = ', dz *Nstep /(V_i_par *Tsim) 
  
  print 'D from', times_r1[0]*dD/r1, ' to ', times_r1[-1]*dD/r1, ' times r1' 
  print 'dD = ', dD, 'm, dD / r1 = ', dD /r1 
  print 'Nzr_lo = ', Nzr_lo, ', Nzr_hi = ', Nzr_hi, ', Nconds_z = ', Nconds_z, '\n'  
  
  t0 = time.time() 
  
  for iD in np.arange(Nsamp_D): 
    if iD%10 == 0: print "iD = ", iD 
    #Nstep = 8000 
    #Nconds_z = ...
    #if iD > 1: 
      #Nstep = 1000 
      #Nconds_z = ...
    
    #z_ini  = np.linspace(Nzr_lo *r_L, Nzr_hi *r_L, Nconds_z) # 
    #F_time_aved = 0.0 *z_ini 
    #dz = z_ini[1] -z_ini[0]
    
    for iz in np.arange(Nconds_z): 
      Fz_arr = magDyn1D_leapfr02(z_ini[iz], -V_i_par, Zrcc, D[iD], Tsim, Nstep)[4] 
      F_free_arr = unpert1D(z_ini[iz], -V_i_par, Zrcc, D[iD], Tsim, Nstep) 
      F_fr = (Fz_arr -F_free_arr) *m_el /1.60218e-19  # eV / m = 1.60218e-19 N  
      F_time_aved[iz] = np.sum(F_fr[1:-1]) /np.float64(Nstep) 
    F_line_integr[iD] = 0.5 *dz *np.sum( F_time_aved[0:Nconds_z-1] +F_time_aved[1:Nconds_z] )
  
  # If necessary, manually overwrite expensive-to-compute elements with pre-computed values: 
  F_line_integr[0] = -5.87602039838814e-05 # computed separately for frac_delta_e_par = .01 
  #F_line_integr[0] = -0.0062005566907 # computed separately for frac_delta_e_par = 1.55 
  F_line_integr[1] = -8.02939724089184e-05 # computed separately for frac_delta_e_par = .01 
  #F_line_integr[2] = -0.00601952418081 # computed separately for frac_delta_e_par = 1.55 
  #F_line_integr[0] = -0.00427158334585 # computed separately for frac_delta_e_par = 2.55 
  t1 = time.time() 
  
  print 'Time in the main loop: ', t1-t0, ' sec.', '\n' 
  
  #print 'Line integral of Fpar = ', 0.5 *dz *np.sum( F_time_aved[0:Nconds_z-1] +F_time_aved[1:Nconds_z] ) , '\n' 
  print 'dD = ', dD, 'm, np.sum(F_line_integr) = ', np.sum(F_line_integr) 
  print 'n_e *2 *pi *np.sum(D*F_line_integr*dD) = ', n_e *2. *np.pi *dD *np.sum(D *F_line_integr), '\n' 
  
  np.savetxt('D_Fintegr_DtimesFintegr.txt', zip(D, F_line_integr, D *F_line_integr), fmt="%0.18e") 
  
  DF = D *F_line_integr 
  r = -(np.log(-DF[-1])  -np.log(-DF[-6])) / (D[-1] -D[-6]) 
  print 'r = ', r 

  S_corrn = n_e *2. *np.pi *dD *DF[-1] / (np.exp(r *dD) -1.) 
  print 'Correction to F = ', S_corrn 
  
  #plt.plot(z_ini, F_time_aved, 'b-') 
  #plt.plot(z_ini, 0.0*z_ini, 'g--') 
  plt.plot(D/r1, F_line_integr, 'bo') 
  plt.plot(D/r1, F_line_integr, 'b-')
  plt.plot(D/r1, 0.0 *D/r1, 'k--') 
  plt.xlim(0.0, ) 
  plt.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
  plt.xlabel('$D /r_1$')
  plt.ylabel('Line integral of $F_{\parallel}$')
  plt.savefig('D_Fintegr.pdf') 
  plt.show() 
  
  plt.plot(D/r1, D *F_line_integr, 'ro') 
  #plt.plot(D, D *F_line_integr, 'r-') 
  plt.plot(D/r1, 0.0 *D/r1, 'g--') 
  plt.xlim(0.0, ) 
  plt.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
  plt.xlabel('$D /r_1$')
  plt.ylabel('D * Line integral of $F_{\parallel}$')
  plt.savefig('D_DtimesFintegr.pdf') 
  plt.show() 
  
  plt.plot(D[-6:D.size], np.log( -DF[-6:D.size] ), 'ro')
  plt.xlabel('D (m)')
  plt.ylabel('$\ln ($|D * Line integral of $F_{\parallel} |)$')
  plt.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
  #plt.savefig('DtimesFintegr_tail_scaling.pdf')
  plt.show() 



def unpert1D(z_ini, v0, C0, D0, T, Nstep): 
  
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  t_arr[-1] = T 
  t_arr[1:Nstep+1] = 0.5*dt +np.linspace(0., T -dt, Nstep) 
  z_free_arr = z_ini +v0 *t_arr 
  r = np.sqrt(D0*D0 +z_free_arr*z_free_arr) 
  F_free_arr = C0 *z_free_arr /(r*r*r)  # from the electron on the ion, hence the "+" sign 
  
  return F_free_arr 




def magDyn1D_leapfr02(z_ini, v0, C0, D0, T, Nstep): 
  dt = T /np.float64(Nstep) 
  t_arr = np.zeros(Nstep+2, dtype=np.float64) 
  z_arr = 1.0*t_arr  
  v_arr = np.zeros(Nstep+1, dtype=np.float64) 
  t_v_arr = np.linspace(0.0, T, Nstep+1) 
  Fz_arr = 1.0*t_arr  # from electron on the ion 
  
  z = z_ini 
  vz = v0 
  #print v0
  
  z_arr[0] = z 
  v_arr[0] = vz 
  r = np.sqrt(D0*D0 +z*z) 
  r3 = r*r*r
  Fz_arr[0] = C0 *z / r3  #  "+", because the force from electron on the ion 
  
  z += vz * dt/2. 
  
  t_arr[1] = dt/2. 
  z_arr[1] = z 
  r = np.sqrt(D0*D0 +z*z) 
  r3 = r*r*r
  Fz_arr[1] = C0 *z / r3
  
  for istep in np.arange(Nstep-1):
    vz += -dt *C0 *z / r3  # C0 > 0 
    z += vz *dt 
    r = np.sqrt(D0*D0 +z*z) 
    r3 = r*r*r 
    
    t_arr[2+istep] = (1.5 +istep) *dt 
    z_arr[2+istep] = z 
    Fz_arr[2+istep] = C0 *z / r3 
    v_arr[1+istep] = vz 
  
  vz += -dt *C0 *z / r3  # C0 > 0 
  z += vz *dt/2. 
  
  t_arr[Nstep+1] = T  
  z_arr[Nstep+1] = z 
  r = np.sqrt(D0*D0 +z*z) 
  r3 = r*r*r
  Fz_arr[Nstep+1] = C0 *z / r3 
  v_arr[Nstep] = vz 
  
  return z, vz, t_arr, z_arr, Fz_arr, t_v_arr, v_arr   



if __name__=="__main__":
  import numpy as np 
  import matplotlib.pyplot as plt 
  import time 
  main() 
  
