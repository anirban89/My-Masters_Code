
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
path2 = home+'/LDinfBeta1test1';
path3 = home+'/LDinfBeta3test1';
path4 = home+'/LDinfBeta5test1';
path1 = home+'/LD1Beta0test1';
path2 = home+'/LD1Beta1test1';
path3 = home+'/LD1Beta3test1';
path4 = home+'/LD1Beta5test1';

#path1 = home+'/PVRampPeriodicLDinfBeta0';
#path2 = home+'/PVRampPeriodicLDinfBeta5';
#path3 = home+'/PVRampPeriodicLD1Beta0';
#path4 = home+'/PVRampPeriodicLD1Beta5';

#path1 = home+'/SymJetsLDinfBeta0';
#path2 = home+'/SymJetsLDinfBeta5';
#path3 = home+'/SymJetsLD1Beta0';
#path4 = home+'/SymJetsLD1Beta5';

#LD = 1;


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
#N = 128;
t1 = 1;
T_tot = 200 ; 
#T_tot = 500 ;
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
tx = [2, 10, 50, 100];
Ex = zeros((4));
Ax = zeros((4));
figure(figsize= (14,6));

for i in arange(N_t) :
	ik = (i+1)*1
	Omega = np.load(path1+'/PV'+str(ik)+'.npy');
	psi_n = InvPotentialVorticity(Omega);
	u_n = -partialY(psi_n);
	v_n = partialX(psi_n);
	Energy[i] = mean(u_n*u_n + v_n*v_n) ;
	Enstrophy[i] = mean(Omega * Omega);
	Anisotropy[i] = mean(u_n*u_n)/(mean(u_n*u_n) + mean(v_n*v_n))
	EnsEnFr[i] = Enstrophy[i]/Energy[i];

Ex = [EnsEnFr[2], EnsEnFr[10], EnsEnFr[50],EnsEnFr[100]]
Ax = [Anisotropy[2], Anisotropy[10], Anisotropy[50], Anisotropy[100]]

subplot(121)

p1 = plot(tn,EnsEnFr,'r')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Dirichlet Quotient, $ \Lambda (t)$ ',fontsize=16)


subplot(122)
p1 = plot(tn,Anisotropy,'r')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Anisotropy, $\mathcal{R}(t)$',fontsize=16)

hold(True)

for i in arange(N_t) :
	ik = (i+1)*1
	Omega = np.load(path2+'/PV'+str(ik)+'.npy');
	psi_n = InvPotentialVorticity(Omega);
	u_n = -partialY(psi_n);
	v_n = partialX(psi_n);
	Energy[i] = mean(u_n*u_n + v_n*v_n) ;
	Enstrophy[i] = mean(Omega * Omega);
	Anisotropy[i] = mean(u_n*u_n)/(mean(u_n*u_n) + mean(v_n*v_n))
	EnsEnFr[i] = Enstrophy[i]/Energy[i];

Ex = [EnsEnFr[2], EnsEnFr[10], EnsEnFr[50],EnsEnFr[100]]
Ax = [Anisotropy[2], Anisotropy[10], Anisotropy[50], Anisotropy[100]]
subplot(121)

p2 = plot(tn,EnsEnFr,'b')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Dirichlet Quotient, $ \Lambda (t)$ ',fontsize=16)


subplot(122)
p2 = plot(tn,Anisotropy,'b')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Anisotropy, $\mathcal{R}(t)$',fontsize=16)
#hold(True)

for i in arange(N_t) :
	ik = (i+1)*1
	Omega = np.load(path3+'/PV'+str(ik)+'.npy');
	psi_n = InvPotentialVorticity(Omega);
	u_n = -partialY(psi_n);
	v_n = partialX(psi_n);
	Energy[i] = mean(u_n*u_n + v_n*v_n) ;
	Enstrophy[i] = mean(Omega * Omega);
	Anisotropy[i] = mean(u_n*u_n)/(mean(u_n*u_n) + mean(v_n*v_n))
	EnsEnFr[i] = Enstrophy[i]/Energy[i];

Ex = [EnsEnFr[2], EnsEnFr[10], EnsEnFr[50],EnsEnFr[100]]
Ax = [Anisotropy[2], Anisotropy[10], Anisotropy[50], Anisotropy[100]]
subplot(121)

p3 = plot(tn,EnsEnFr,'g')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Dirichlet Quotient, $ \Lambda (t)$ ',fontsize=16)


subplot(122)
p3 = plot(tn,Anisotropy,'g')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Anisotropy, $\mathcal{R}(t)$',fontsize=16)

for i in arange(N_t) :
	ik = (i+1)*1
	Omega = np.load(path4+'/PV'+str(ik)+'.npy');
	psi_n = InvPotentialVorticity(Omega);
	u_n = -partialY(psi_n);
	v_n = partialX(psi_n);
	Energy[i] = mean(u_n*u_n + v_n*v_n) ;
	Enstrophy[i] = mean(Omega * Omega);
	Anisotropy[i] = mean(u_n*u_n)/(mean(u_n*u_n) + mean(v_n*v_n))
	EnsEnFr[i] = Enstrophy[i]/Energy[i];

Ex = [EnsEnFr[2], EnsEnFr[10], EnsEnFr[50],EnsEnFr[100]]
Ax = [Anisotropy[2], Anisotropy[10], Anisotropy[50], Anisotropy[100]]
subplot(121)

p4 = plot(tn,EnsEnFr,'k')
xlabel(r'Time, $t$',fontsize=16)
ylabel(r'Dirichlet Quotient, $ \Lambda (t)$ ',fontsize=16)
legend([p4, p3, p2, p1], [r'$\beta = 5$', r'$\beta = 3$', r'$\beta = 1$', r'$\beta = 0$'])
#legend([p2, p1, p4, p3], [r'$L_D = \infty$, $\beta = 5$', r'$L_D = \infty$, $\beta = 0$', r'$L_D = 1$, $\beta = 5$', r'$L_D = 1$, $\beta = 0$'])
#legend([p2, p4, p3], [r'$L_D = \infty$, $\beta = 5$', r'$L_D = 1$, $\beta = 5$', r'$L_D = 1$, $\beta = 0$'])
subplot(122)
p4 = plot(tn,Anisotropy,'k')
xlabel(r'Time, $t$',fontsize=16)
ylabel( r'Anisotropy, $\mathcal{R}(t)$',fontsize=16)
#legend([p4, p3, p2, p1], [r'$\beta = 5$', r'$\beta = 3$', r'$\beta = 1$', r'$\beta = 0$'])
#legend([p2, p1, p4, p3], [r'$L_D = \infty$, $\beta = 5$', r'$L_D = \infty$, $\beta = 0$', r'$L_D = 1$, $\beta = 5$', r'$L_D = 1$, $\beta = 0$'])



savefig('/home/anirban/LD_1_allBetaRhIC_loff.png');
#savefig('/home/anirban/SymJets2_loff.png');
#savefig('/home/joymm/anirban/Bickley-Vars_LD'+str(LD)+'Beta-'+str(Beta)+'.png');

show()
