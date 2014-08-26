from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;
import sys;
import pylab;
from os.path import expanduser;

home = expanduser("~");

if(len(sys.argv) > 1):
    
        alphaVapor = float(sys.argv[1]);
        betaVapor = float(sys.argv[2]);
        print "alphaQ = ", str(alphaVapor);
        print "betaQ = ", str(betaVapor);
else:
    print 'Insufficient arguments';
    exit();

N = 256;
h = 2*pi/N;

tn = 0;
sec = 1;

initPath = home+'//TwoLayerLinear';

if alphaVapor<0:
    initPath = home+'//TwoLayerNegative';

path = initPath+str(int(10*abs(alphaVapor)))+str(int(10*betaVapor));
figPath = home+'/'+str(int(10*abs(alphaVapor)))+str(int(10*betaVapor));

dx = 2*pi/(N-1);
dy = 2*pi/(N-1);
qs = zeros((N,N),dtype=double);
q0 = 7;

betaMax = q0/(2*pi);
alphaMax = q0/(2*pi);

beta = betaVapor*betaMax;
alpha = (alphaVapor*alphaMax);


#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
#initPath = '/Data/';
vorticityPath1 = '/eta1_data_'
vorticityPath2 = '/eta2_data_'
vaporPath = '/vapor_data_'

k = arange(0,N/3);
def deriveVelocities(omega1, omega2):
    '''
    omega1 = laplacian(psi1) + (psi2-psi1)
    omega2 = laplacian(psi2) + (psi1-psi2)
    '''
    sumPsi = InvLaplacian(omega1+omega2);
    diffPsi = InvPotentialVorticityTwoLayer(omega2 - omega1);

    psi1 = (sumPsi + diffPsi)/2.0;
    psi2 = (sumPsi - diffPsi)/2.0;

    return(psi1,psi2);


def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        val = fft(field[:,i]);
        power = power + abs(val*val.conj());


    power = power/y;
    return log(power[0:N/3]);

def calcVaporPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        val = fft(field[:,i]);
        power = power + abs(val*val.conj());

    power = power/y;
    return log(power[1:N/2]);



hold(True);

Xrange=arange(0,2*pi,2*pi/256)
Yrange=arange(0,2*pi,2*pi/256)

data1 =  np.load(path+'/vapor_data_1.npy');
data1[data1<0] = 0;
CS = contour(Xrange,Yrange,(data1),colors='k')
#imshow((data),extent=(0,2*pi,0,2*pi))
#colorbar();
data5 =  np.load(path+'/vapor_data_50.npy');
data5[data5<0] = 0;
#imshow((data5-data1),extent=(0,2*pi,0,2*pi),vmin=-0.2,vmax=0.3)
#colorbar();
CS = contour(Xrange,Yrange,(data5),colors='b')
data10 =  np.load(path+'/vapor_data_200.npy');
data10[data10<0] = 0;
#imshow((data10-data1),extent=(0,2*pi,0,2*pi),vmin=-0.2,vmax=0.3)
#colorbar();
contour(Xrange,Yrange,(data10),colors='r')

xlabel(r' $ \leftarrow$ West ',fontsize=16)

ylabel(r"South $\leftarrow$ ",fontsize=16)


savefig('/home/joymm/lowResRuns/contours-'+str(int(10*alphaVapor))+str(int(10*betaVapor))+'.jpg',dpi=400);

show();


