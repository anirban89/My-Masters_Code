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
path = home+'/PVStepPeriodicLDinfBeta5';

if(len(sys.argv) > 1):
    if(len(sys.argv) == 2):
    
        print "Resuming from";
        resume = int(sys.argv[1]);
        print resume;

Beta = 5;
N = 256;
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

#for i in arange(N):
#  
#       if (i<N/4):
#           omega[:,i]=0.5;
#       elif (i >= N/4 and i<N/2):
#           omega[:,i]=1.5;
#       elif (i >= N/2 and i<3*N/4):
#           omega[:,i]=2.5;
#       else:
#           omega[:,i]=3.5;

"""
PERTURBATION
"""


for i in arange(N):
    for j in arange(N):
        psi_p[i,j] = cos(x[i]+0.3) + 0.9*sin(3*(y[j]+1.8) + 2*x[i]) + 0.87*sin(4*(x[i]-0.7) + (y[j]+0.4) + 0.815) + 0.8*sin(5*(x[i]-4.3)+0.333) + 0.7*sin(7*y[j]+0.111)
omega1_p = partialX(psi_p,2) + partialY(psi_p,2) ;
for j in arange(N):
	omega_p[:,j] = omega1_p[:,j]
"""
BASE FLOW
"""
n = N;
for i in arange(n):
  
       if (i<n/6):
           omega[:,i]=0;
       elif (i >= n/6 and i<n/3):
           omega[:,i]=3;
       elif (i >= n/3 and i<2*n/3):
           omega[:,i]=6;
       elif (i >= 2*n/3 and i<5*n/6):
           omega[:,i]=3;	
       else:
           omega[:,i]=0;

omega = (omega + 0.005*omega_p);

    

  
"""
for i in arange(N):
    for j in arange(N):
        si[i,j] = cos(x[i]+0.3) + 0.9*sin(3*(y[j]+1.8) + 2*x[i]) + 0.87*sin(4*(x[i]-0.7) + (y[j]+0.4) + 0.815) + 0.8*sin(5*(x[i]-4.3)+0.333) + 0.7*sin(7*y[j]+0.111);
	
omega1 = partialX(partialX(si)) + partialY(partialY(si)) - si ;
for j in arange(N):
	omega[:,j] = omega1[:,j] + Beta*y[j];
"""

def domegadt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v, forcing) = args;

    return real(-(u*partialX(omega) + v*partialY(omega) + Beta*v + forcing));


def calcMax_omega(fields, args):
    omegaNew = fields;

    psi = InvPotentialVorticity(omegaNew);
    u = -partialY(psi);
    v = partialX(psi);
	
    maximum = \
        sqrt(amax((amax(u*u),amax(v*v))));
    #print 'Max vel: ', maximum;
    return maximum;

def dSystemdt(dt,fields,args):
    omegaNew = fields;

    forcing = args;
    psi = InvPotentialVorticity(omegaNew);
    u = -partialY(psi);
    v = partialX(psi);
    dOmega = domegadt(dt, omegaNew, (u,v,forcing));
    return(dOmega);


def diffusion(dt, fields, args):
    omegaNew = fields;

    diffOmega = spectralDiffusion(dt, omegaNew, args);

    return(diffOmega);


stepper = integrator(h, [dSystemdt], \
                [diffusion], calcMax_omega,1);

ion();
omega_n = omega;
psi_n = InvPotentialVorticity(omega_n);
u_n = -partialY(psi_n);
v_n = partialX(psi_n);
Energy = mean(u_n*u_n) + mean(v_n*v_n) ;
Enstrophy = mean(omega_n * omega_n);
Ttt = 0.0;
en = 0.0;
#figure(figsize= (10,8));
figure(figsize= (16,5));
subplot(121)

prev_t = 0;
tn2 = 0;
sec = 0;
#jet();

if(resume > 0):
   #print 'Resuming from step: ',resume;
   omega_n = \
       np.load(path+'/eta_data_'+str(resume)+'.npy').transpose();
   sec = resume+10;
   prev_t = resume;



handle1 = imshow(omega_n.transpose(), vmin=-15,vmax=1, aspect='auto');
c4 = colorbar();
draw();

#figure(figsize= (10,8));
subplot(122)
handle2 = imshow(psi_n.transpose(), vmin=-15,vmax=1, aspect='auto');
c4 = colorbar();
draw();

#subplot(212)
#handle3 = subplot(212)

global_time = 0;
forcing = 0.0*rand(N,N);
forcingTime = 1;
tn = 0;
Ttot = 600;
dtx = 1;

Nt = Ttot/dtx;
#Ttt = zeros(Nt);
#en = zeros(Nt);
#for ix in arange(Nt):
en = []
while (sec<Ttot):
    tx = amax(sec);
    en.append(tx);
    #forcing = rand(N,N);
    (dt,fields) = stepper.integrate(tn,omega_n,forcing);
    omega_n = fields;
    psi_n = InvPotentialVorticity(omega_n);
    u_n = -partialY(psi_n);
    v_n = partialX(psi_n);
    Energy = mean(u_n*u_n) + mean(v_n*v_n) ;
    Enstrophy = mean(omega_n * omega_n);
    ens_ene_frac = Enstrophy/Energy ;
    #en[sec] = Energy;
    tn = dt + prev_t;
    prev_t = tn;
    #print '============================='
    #print 'Time Elapsed: ', tn;
    #print '============================='
    print '============================='
    print 'Energy of the system',Energy;
    print '============================='
    print '============================='
    print 'Enstrophy of the system',Enstrophy;
    print '============================='
    print '============================='
    print 'Enstrophy/Energy of the system',ens_ene_frac;
    print '============================='
    if tn>forcingTime:
	#print 'changing forcing';
	forcingTime += 1;
	forcing = 0.0*rand(N,N);
    if(tn > sec):
	#print '==================================================================='
	#print 'Time elapsed: ',tn;
	#print '==================================================================='

	handle1.set_array(omega_n.transpose());
	handle1.autoscale();

	handle2.set_array(psi_n.transpose());
	handle2.autoscale();
	sec = sec+dtx;
	
	#while (sec>1):
	#print  sec , en
	#handle3.plot(en);
	hold();
	draw();
	
	np.save(path+'/PV'+str(sec)+'.npy',omega_n);
	savefig(path+'/fig'+str(sec)+'.png');
	
ioff();
