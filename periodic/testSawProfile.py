from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
import scipy as sp;
import pylab as m;
import sys;
import numpy as np;
from os.path import expanduser

N = 512;
h = 2*pi/N;

home = expanduser("~");


path = home+'/PVRamptestInitialization';


Beta = 0;
x = linspace(0,2*pi, N+1);
x = x[0:-1];
y = linspace(0,2*pi,N+1);
y = y[0:-1];
y1 = linspace(0,2*pi,N+1);
y1 = y1[0:-1];
Uimp = -0.0;

omega = 0.5*rand(N,N);
omega_p = 0.5*rand(N,N);
si = zeros((N,N));
for i in arange(N):
    for j in arange(N):
        si[i,j] = cos(x[i]+0.3) + 0.9*sin(3*(y1[j]+1.8) + 2*x[i]) + 0.87*sin(4*(x[i]-0.7) + (y1[j]+0.4) + 0.815) + 0.8*sin(5*(x[i]-4.3)+0.333) + 0.7*sin(7*y1[j]+0.111);



n = N ;

"""
for i in arange(n):
  
       if (i<n/6):
           omega[:,i]=0;
       elif (i >= n/6 and i<n/3):
           omega[:,i]=-10*(y[i] - pi/6);
       elif (i >= n/3 and i<n/2):
           omega[:,i]=-10*(y[i] - pi/3);
       elif (i >= n/2 and i<2*n/3):
           omega[:,i]=-10*(y[i] - pi/2);
       elif (i >= 2*n/3 and i<5*n/6):
           omega[:,i]=-10*(y[i] - 2*pi/3);	
       else:
           omega[:,i]=0;
"""

for i in arange(n):
  
       if (i<n/6):
           omega[:,i]=-50*(y[i]);
       elif (i >= n/6 and i<n/4):
           omega[:,i]=100*(y[i] - pi/2);
       elif (i >= n/4 and i<5*n/12):
           omega[:,i]=-50*(y[i] - pi/2);
       elif (i >= 5*n/12 and i<n/2):
           omega[:,i]=100*(y[i] - pi);
       elif (i >= n/2 and i<2*n/3):
           omega[:,i]=-50*(y[i] - pi);
       elif (i >= 2*n/3 and i<3*n/4):
           omega[:,i]=100*(y[i] - 3*pi/2);
       elif (i >= 3*n/4 and i<11*n/12):
           omega[:,i]=-50*(y[i] - 3*pi/2);	
       else:
           omega[:,i]=100*(y[i] - 2*pi);


omega = (omega + 0.001*omega_p*exp(-pow(y-pi/2,2)));
si = InvPotentialVorticity(omega);
omega_p = partialX(si,2) + partialY(si,2);
#figure( figsize = (18,5)
subplot(131)
imshow((omega).transpose())
colorbar()
subplot(132)
imshow((omega_p).transpose())
colorbar();
subplot(133)
imshow((si).transpose())
colorbar();
print amax(omega - omega_p);

np.save(path+'/SawPVInitialPeriodic2pi.npy',si);
savefig(path+'/figInitializeSaw.png');

show();



