from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 256;
kappa = 0.2;
initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
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
corr = zeros((N));

spec = zeros((510,N));

for i in range(0,510):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();
    psi = InvPotentialVorticity(pvData);
    u = -partialY(psi)
    print shape(pvData);
    print 'reading',vorticityPath+str(i)+'.npy';
    print shape(average(u,-1));
    spec[i,:] = fft(average(u,-1));

spacetime = zeros((510,N));

for i in arange(N):
    spacetime[:,i] = fft(spec[:,i]);


imshow(abs(spacetime*conj(spacetime))[0:205,0:N/2]);
show()

