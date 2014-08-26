from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
from math import *
import scipy as sp;
import pylab as m;
import sys;
import numpy as np;

Beta = 0;
N = 128;
h = 2*pi/N;
x = linspace(0,2*pi, N+1);
x = x[0:-1];
y = linspace(0,2*pi, N+1);
y = y[0:-1];
y1 = linspace(-5,5, N+1);
y1 = y1[0:-1]; 
Uimp = -0.0;
omega1 = zeros((N,N));
omega1_p = zeros((N,N));
omega = 0.5*rand(N,N);
si = 0.5*rand(N,N);
omega_p = rand(N,N);
psi_p = rand(N,N);

for i in arange(N):
    #si[:,i] = -tanh(y1[i]);
    #omega1[:,i] = (3-2*pow((tanh(y1[i])),2))*tanh(y1[i]);
    #omega1[:,i] = si[:,i];
    si[:,i] = sin(y[i]);
    si[:,i] = y[i];

#omega =  partialY(si,2) + partialX(si,2) - si ;
omega = partialY(si);

imshow(omega.transpose())
colorbar()
show()
