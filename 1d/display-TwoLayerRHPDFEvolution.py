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
    if(len(sys.argv) == 4):
    
        alphaVapor = float(sys.argv[1]);
        betaVapor = float(sys.argv[2]);
        rangeVals = int(sys.argv[3]);
        print "alphaQ = ", str(alphaVapor);
        print "betaQ = ", str(betaVapor);
else:
    print 'Insufficient arguments';
    exit();

N = 256;
h = 2*pi/N;

tn = 0;
sec = 1;

initPath = home+'/TwoLayer';

if alphaVapor<0:
    initPath = home+'/TwoLayerNegative';

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

if(alphaVapor > 0):
    for i in arange(N):
        for j in arange(N):
            qs[N-i-1,j] = q0 - alpha*dx*(i) - beta*dy*j;
elif(alphaVapor < 0):
    for i in arange(N):
        for j in arange(N):
            qs[i,j] = q0 - alpha*dx*(i) - beta*dy*j;
else:
    for i in arange(N):
        for j in arange(N):
            qs[i,j] = q0 - beta*dy*j;



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


# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
wavenumber = arange(0,N/3);
wavenumberVapor = arange(0,N/2);

vaporHovmoller = zeros((rangeVals,N));
pvHovmoller = zeros((rangeVals,N));
condensation = zeros(rangeVals);
vaporVals = zeros(rangeVals);

#ion();
hold(True);

#axes.set_autoscale_on(True)
#axes.autoscale_view(True,True,True)

bins = linspace(0,1,10);
centre = (bins[1:len(bins)] + bins[0:len(bins)-1])/2.0;
handle, = plot(centre,zeros(len(centre)));

maxVal = 0;
for i in arange(10,rangeVals,1000):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
#    pvData1 = np.load(path+vorticityPath1+str(i)+'.npy').transpose();
#    pvData2 = np.load(path+vorticityPath2+str(i)+'.npy').transpose();
    average = zeros(len(centre));
    for j in arange(i, i+1000, 10):
        print 'reading time step:', str(j);
        vaporData = np.load(path+vaporPath+str(j)+'.npy').transpose();
        rh = (vaporData+qs)/qs;
        rh = rh[rh<=1];
        hist, edges = histogram(rh, bins, density=True);
        print len(hist);
        print len(edges);
        print len(centre);
        average += hist;

    #handle.set_data(centre,average/100.0);
    plot(centre,average/100.0);
    if(maxVal < max(average)):
        maxVal = max(average);
#    axes.relim();
#    axes.autoscale_view(True,True,True);
#    draw();


ylim([0, maxVal/100.0]);
xlabel('RH');
ylabel('PDF');
title('Evolution of PDF over time for alpha = '+str(alphaVapor));
show();
#    vaporData[vaporData < 0] = 0;
#    vaporData[i/10,:,:] = (qs+vaporData[i/10,:,:]);
#    actual = qs;
#    rh = actual/qs;
#    print amax(rh);





