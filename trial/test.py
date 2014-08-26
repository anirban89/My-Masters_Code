from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
import scipy as sp;
import pylab as m;
import sys;
import numpy as np;
#from os.path import expanduser
#home = expanduser("~");
N = 128;
h = 2*pi/N;

u = 0.5* rand(N,N);
v = 0.5* rand(N,N);


for i in arange(N):
    for j in arange(N):
        u[i,j] = 2*cos(i*h)*sin(3*j*h);
	v[i,j] = 4*sin(i*h)*cos(3*j*h);

handle = imshow(u.transpose())

figure()
imshow(v.transpose())

#figure(figsize= (10,8));
#handle = imshow(u.transpose(), vmin=-2, vmax=0);
#draw();

meanu = mean(u);
u_rms = sqrt(mean(u*u));
v_rms = sqrt(mean(v*v));
Energy = u_rms*u_rms + v_rms*v_rms ;
print u_rms, v_rms , Energy 

show()

