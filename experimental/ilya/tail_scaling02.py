import numpy as np
import matplotlib.pyplot as plt 

D, F, DF = np.loadtxt("D_Fintegr_DtimesFintegr.txt", skiprows=0 , unpack=True)

n_e = 2.0e+15 
dD = D[2]-D[1] #4.923725109213487e-07 # r1  

#print n_e *2. *np.pi *dD *np.sum(DF) 

plt.plot(D[-6:D.size], np.log( -DF[-6:D.size] ), 'ro')
plt.xlabel('D (m)')
plt.ylabel('$\ln ($|D * Line integral of $F_{\parallel} |)$')
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
#plt.savefig('DtimesFintegr_tail_scaling.pdf')
plt.show() 

#r = -(np.log(-DF[-1])  -np.log(-DF[-10])) / (D[-1] -D[-10]) 
r = -(np.log(-DF[-1])  -np.log(-DF[-6])) / (D[-1] -D[-6])
print 'r = ', r 

S_corrn = n_e *2. *np.pi *dD *DF[-1] / (np.exp(r *dD) -1.) 
print 'Correction to F = ', S_corrn

