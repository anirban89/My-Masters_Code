from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 500;
kappa = 0.4;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/SlowCondensation/';
vorticityPath = 'Vapor'+str(kappa)+'/eta_data_'
vaporPath = 'Vapor'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
forcing = zeros((N))
const = zeros((N));
vel = zeros((N));
water = zeros((N));
for i in range(N):

    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();
    #vaporData = np.load(initPath+vaporPath+str(i)+'.npy').transpose();
    print 'reading',vorticityPath+str(i)+'.npy';
    psi = InvPotentialVorticity(pvData);
    v = partialX(psi);

    vel[i] = average(abs(v));


avg = zeros(480);
j = 0;
for i in range(480):
    avg[j] = log(average(vel[i:i+5])+0.01);

plot(vel, label='V velocity');
legend(loc='lower left');
title('Value of V over time for k = '+str(kappa));
show();
