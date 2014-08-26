from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
from math import *
import scipy as sp;
import pylab as m;
import sys;
import numpy as np;
from os.path import expanduser

home = expanduser("~");

resume = 0;
path = home+'/LDinfBeta0test1';

N = 128;
t1 = 2;
T_tot = 500 ; 
N_t = 500 ;
tn = linspace(t1,T_tot, N_t);
Omega = zeros((N,N));
psi_n = zeros((N,N));
Energy = zeros((N_t));
Enstrophy = zeros((N_t));
EnsEnFr = zeros((N_t));
for i in arange(N_t) :
	ik = (i+1)*1
	Omega = np.load(path+'/PV'+str(ik)+'.npy');
	psi_n = InvPotentialVorticity(Omega);
	u_n = -partialY(psi_n);
	v_n = partialX(psi_n);
	Energy[i] = mean(u_n*u_n) + mean(v_n*v_n) ;
	Enstrophy[i] = mean(Omega * Omega);

EnsEnFr = Enstrophy/Energy;
figure(figsize= (16,5));
subplot(131)
plot(tn,Energy)
subplot(132)
plot(tn,Enstrophy)
subplot(133)
plot(tn,EnsEnFr)
show()
