from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
import scipy as sp;
import pylab as m;
import sys;
import numpy as np;
from os.path import expanduser

home = expanduser("~");

resume = 0;
path = home+'/Betaplanechannel1';

if(len(sys.argv) > 1):
    if(len(sys.argv) == 2):
    
        print "Resuming from";
        resume = int(sys.argv[1]);
        print resume;

N = 256;
h = 2*pi/N;


x = linspace(0,2*pi, N+1);
x = x[0:-1];
y = linspace(0,pi,N/2+1);
y = y[1:-1];

Uimp = -0.0;

omega = 0.5*rand(N,N/2-1);
for i in arange(N):
    omega[i,:] = cos(2*x[i])*sin(2*y)*exp(-pow(y-pi/2,2));

#for i in arange(N):
#    for j in arange(N):
#        omega[i,j] = 2*cos(i*h)*sin(3*j*h);
#

def domegadt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v, forcing) = args;

    return real(-(u*partialX(omega) + v*partialYSine(omega) + v + forcing));


def calcMax_omega(fields, args):
    omegaNew = fields;

    psi = InvPotentialVorticitySine(omegaNew);
    u = -partialYSine(psi);
    v = partialX(psi);
	
    maximum = \
        sqrt(amax((amax(u*u),amax(v*v))));
    print 'Max vel: ', maximum;
    return maximum;

def dSystemdt(dt,fields,args):
    omegaNew = fields;

    forcing = args;
    psi = InvPotentialVorticitySine(omegaNew);
    u = -partialYSine(psi);
    v = partialX(psi);
    dOmega = domegadt(dt, omegaNew, (u,v,forcing));
    #return(zeros((N,N/2-1)));
    return(dOmega);


def diffusion(dt, fields, args):

    return fields;
    leny = fields.shape[1];    
    lenx = fields.shape[0];
    newLen = 2*(leny+1);
    omegaNew = zeros((lenx,newLen));
    for i in arange(leny):
        omegaNew[:,i+1] = fields[:,i];
        omegaNew[:,newLen-i-1] = -fields[:,i];
    
    diffOmega = spectralDiffusion(dt, omegaNew, args);

    return(diffOmega[:,0:newLen/2-1]);


stepper = integrator(h, [dSystemdt], \
                [diffusion], calcMax_omega,1);

ion();
omega_n = omega;

figure(figsize= (10,8));


prev_t = 0;
tn2 = 0;
sec = 0;
#jet();

if(resume > 0):
   print 'Resuming from step: ',resume;
   omega_n = \
       np.load(path+'/eta_data_'+str(resume)+'.npy').transpose();
   sec = resume+10;
   prev_t = resume;



handle1 = imshow(omega_n.transpose(), vmin=-15,vmax=1, aspect='auto');
c4 = colorbar();
draw();

global_time = 0;
forcing = 0.5*rand(N,N/2 - 1);
forcingTime = 1;
tn = 0;

while (sec<100000):

    #forcing = rand(N,N);
    (dt,fields) = stepper.integrate(tn,omega_n,forcing);

    omega_n = fields;

    tn = dt + prev_t;
    prev_t = tn;
    print '============================='
    print 'Time Elapsed: ', tn;
    print '============================='

    if tn>forcingTime:
        print 'changing forcing';
        forcingTime += 1;
        forcing = 0.5*rand(N,N/2 - 1);

    if(tn > sec):


        print '==================================================================='

        print 'Time elapsed: ',tn;

        print '==================================================================='

        handle1.set_array(omega_n.transpose());
        handle1.autoscale();
        sec = sec+.1;
        draw();
#	np.save(path+'/PV'+str(sec)+'.npy',omega_n);
        #savefig(path+'/fig'+str(sec)+'.png');

ioff();

