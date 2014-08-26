
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
from scipy import *

home = expanduser("~");

resume = 0;

path1 = home+'/LDinfBeta0test1';
#path2 = home+'/LDinfBeta1test1';
#path3 = home+'/LDinfBeta3test1';
#path4 = home+'/LDinfBeta5test1';
#path1 = home+'/LD1Beta0test1';
#path2 = home+'/LD1Beta1test1';
#path3 = home+'/LD1Beta3test1';
#path4 = home+'/LD1Beta5test1';
#path1 = home+'/PVSawPeriodicLDinfBeta0';
#path2 = home+'/PVSawPeriodicLDinfBeta5';
#path3 = home+'/PVSawPeriodicLD1Beta0';
#path4 = home+'/PVSawPeriodicLD1Beta5';
#LD = 1;


Beta = 0;

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        val = fft(field[:,i]);
        power = power + abs(val*val.conj());


    power = power/y;
    print log(power)
    return log(power[0:N/3]);



N = 256;
t1 = 1;
T_tot = 200 ; 
N_t = T_tot ;
tn = linspace(t1,T_tot, N_t);
Omega = zeros((N,N));
psi_n = zeros((N,N));
Energy = zeros((N_t));
Enstrophy = zeros((N_t));
EnsEnFr = zeros((N_t));
Anisotropy = zeros((N_t));
U = zeros((N,N));
Uyy = zeros((N,N));
RK = zeros((N,N));


tx = 100;

Omega = np.load(path1+'/PV'+str(tx)+'.npy');
psi_n = InvPotentialVorticity(Omega);
psi_nf = flipud(psi_n.transpose())
U = -partialY(psi_n);
Uf = flipud(U.transpose());
Uyy = partialY(U,2);

RK = Uyy - Beta;
RKf = flipud(RK.transpose())

figure(figsize= (9,6))
imshow(psi_nf,aspect='auto',extent=(0,2*pi,0,2*pi));
title(r' Rayleigh-Kuo criterion, $U'' -\beta $  at t ='+str(tx),fontsize=16)
xlabel(' x ',fontsize=16)

ylabel(' y ',fontsize=16)

colorbar()
savefig('/home/anirban/RhIC-RK-LDinf-Beta-'+str(Beta)+'-t-'+str(tx)+'.png');
