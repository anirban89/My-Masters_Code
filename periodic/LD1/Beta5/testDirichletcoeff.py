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
path = home+'/LD1Beta5test1';

N = 128;
t1 = 2;
T_tot = 100 ; 
N_t = 100 ;
tn = linspace(t1,T_tot, N_t);
Omega = zeros((N,N));
psi_n = zeros((N,N));
Energy = zeros((N_t));
Enstrophy = zeros((N_t));
EnsEnFr = zeros((N_t));
Anisotropy = zeros((N_t));
for i in arange(N_t) :
	ik = (i+1)*1
	Omega = np.load(path+'/PV'+str(ik)+'.npy');
	psi_n = InvPotentialVorticity(Omega);
	u_n = -partialY(psi_n);
	v_n = partialX(psi_n);
	Energy[i] = mean(u_n*u_n + v_n*v_n) ;
	Enstrophy[i] = mean(Omega * Omega);
	Anisotropy[i] = mean(u_n*u_n)/(mean(u_n*u_n) + mean(v_n*v_n))
EnsEnFr = Enstrophy/Energy;
figure(figsize= (10,6));
subplot(221)
plot(tn,Energy)
title(" Energy ",fontsize=16)
subplot(222)
plot(tn,Enstrophy)
title(" Enstrophy ",fontsize=16)
subplot(223)
plot(tn,EnsEnFr)
title(" Enstrophy/Energy ",fontsize=16)
subplot(224)
plot(tn,Anisotropy)
title(" Anisotropy ",fontsize=16)


show();
