#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

LL= 1

dt=np.loadtxt( 'main.log', skiprows=1 )

T= 2*np.pi*0.2

t=dt[:,1] / T;

Ek=dt[:,3] #;  L2=dt[:,4]
E0 = Ek[0]

DEk = Ek/Ek[0] - 1

#plt.axis('scaled')

plt.plot( t , DEk)
#plt.plot( t , np.log10(L2))

plt.xlabel(r'$t/T$')
plt.ylabel(r'$\Delta E_\mathrm{k} /E_0$')

#plt.ylim([-0.4 , 1 ])

plt.savefig( 'Ek.png' )
