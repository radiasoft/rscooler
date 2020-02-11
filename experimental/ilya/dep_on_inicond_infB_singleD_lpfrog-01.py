def main(): 
  
  Tsim = 0.5* 4.0e-10  # s, in the beam frame 
  Nstep = 1000/2  # number of timesteps 
  
  Nconds_z = 12001  #1501 
  
  Z = 1 
  #Bz = 5.0  # T (working in MKS units)
  n_e = 2.0e+15  # m^-3, in the beam frame 
  Delta_e_par = 1.0e+5  # m/s 
  #Delta_e_tr = 4.2e+5  # m/s 
  
  r1 = np.power(3./(4.*np.pi*n_e), 1./3.)  # expectation value of the distance to the nearest electron, assuming isotropic const density 
  
  frac_delta_e_par = .01 #3.05 #0.75  # ratio of the ion velocity to the longit. thermal electron velocity (abs value) 
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
  
  times_r1 = 0.15 #0.145497438 
  D = times_r1 *r1 
  
  T_lin = 2. *np.pi *np.sqrt( D**3 /Zrcc ) 
  
  F_time_aved = 0.0 *z_ini 
  
  print ' ' 
  print "Z = ", Z 
  print 'n_e (beam fr.) = ', n_e, 'm^-3' 
  print "r1 == r_nearest (in 3D) = ", r1, "m" 
  print 'V_i_par = ', V_i_par, 'm/s, Tsim = ', Tsim, 's'  
  print 'V_i_par *Tsim = ', V_i_par *Tsim, 'm,  V_i_par *Tsim /r1 = ', V_i_par *Tsim /r1 
  print 'Tmin = ', T_lin, 's,  0.5*Tmin /Tsim = ', 0.5 *T_lin /Tsim 
  print 'D = ', D, 'm,  D /r1 = ', D /r1  
  print 'times_r1_z_lo = ', Nzr_lo, ', times_r1_z_hi = ', Nzr_hi, ', Nconds_z = ', Nconds_z 
  print 'Nstep = ', Nstep, ', dt = ', Tsim /np.float64(Nstep), 's'
  print 'V_i_par *dt = ', V_i_par *Tsim /np.float64(Nstep), 'm,  dz / (V_i_par *dt)  = ', dz *Nstep /(V_i_par *Tsim), '\n'
  
  t0 = time.time() 
  for iz in np.arange(Nconds_z): 
    Fz_arr = magDyn1D_leapfr02(z_ini[iz], -V_i_par, Zrcc, D, Tsim, Nstep)[4] 
    F_free_arr = unpert1D(z_ini[iz], -V_i_par, Zrcc, D, Tsim, Nstep) 
    F_fr = (Fz_arr -F_free_arr) *m_el /1.60218e-19  # eV / m = 1.60218e-19 N 
    F_time_aved[iz] = np.sum(F_fr[1:-1]) /np.float64(Nstep) 
    
  t1 = time.time() 
  print 'Time in the loop = ', t1 - t0, 's' 
  print 'Line integral of Fpar = ', 0.5 *dz *np.sum( F_time_aved[0:Nconds_z-1] +F_time_aved[1:Nconds_z] ) , '\n' 
  
  plt.plot(z_ini/r1, F_time_aved, 'b-') 
  plt.plot(z_ini/r1, 0.0*z_ini/r1, 'k--') 
  plt.ticklabel_format(axis='both', style='sci', scilimits=(-3,3)) 
  plt.xlim(z_ini[0]/r1, z_ini[-1]/r1) 
  plt.xlabel('$z / r_1$')
  plt.ylabel('Time-averaged $F_{\parallel}(z; D)$')
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
  import time 
  import matplotlib.pyplot as plt 
  main() 
  
