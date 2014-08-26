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
y = linspace(0,pi,N/2+1);
y = y[1:-1];
y1 = linspace(0,2*pi,N/2+2);
y1 = y1[0:-1];
Uimp = -0.0;



Psi = 0.5*rand(N,N)
omega = 0.5*rand(N,N/2-1);
omega_p = 0.5*rand(N,N/2-1);
psi = zeros((N,N/2-1));
si = zeros((N,N/2-1));
for i in arange(N):
    for j in arange(N/2-1):
        si[i,j] = cos(x[i]+0.3) + 0.9*sin(3*(y1[j]+1.8) + 2*x[i]) + 0.87*sin(4*(x[i]-0.7) + (y1[j]+0.4) + 0.815) + 0.8*sin(5*(x[i]-4.3)+0.333) + 0.7*sin(7*y1[j]+0.111);



n = N/2 -1;
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

for i in arange(n):
  
       if (i<n/6):
           psi[:,i]=-5*(y[i]);
       elif (i >= n/6 and i<n/4):
           psi[:,i]=10*(y[i] - pi/4);
       elif (i >= n/4 and i<5*n/12):
           psi[:,i]=-5*(y[i] - pi/4);
       elif (i >= 5*n/12 and i<n/2):
           psi[:,i]=10*(y[i] - pi/2);
       elif (i >= n/2 and i<2*n/3):
           psi[:,i]=-5*(y[i] - pi/2);
       elif (i >= 2*n/3 and i<3*n/4):
           psi[:,i]=10*(y[i] - 3*pi/4);
       elif (i >= 3*n/4 and i<11*n/12):
           psi[:,i]=-5*(y[i] - 3*pi/4);	
       else:
           psi[:,i]=10*(y[i] - pi);
"""
Psi = np.load(path+'/SawPVInitialPeriodic2pi.npy');

psi = Psi[:,0:N/2-1]
si = psi - amin(psi);
omega = partialX(psi,2) + partialYSine(psi,2) ;

#omega = (omega + 0.001*omega_p*exp(-pow(y-pi/2,2)));

#si = InvPotentialVorticitySine(omega);

omega_p = partialX(si,2) + partialYSine(si,2);

subplot(221)
imshow((psi).transpose())
colorbar()

subplot(222)
imshow((omega).transpose())
colorbar()

subplot(223)
imshow((si).transpose())
colorbar();

subplot(224)
imshow((omega_p).transpose())
colorbar();
figure()
plot(y,omega[14,:])

print amax(omega - omega_p);

show();



