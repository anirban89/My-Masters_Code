from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
import scipy as sp;
import pylab as m;
import sys;
import numpy as np;

N = 128;

x = arange(100);
z = zeros((10,10))
h = 2*pi/100.0;
y = sin(x*h);

fig = plot(x,y)

print z
show();
