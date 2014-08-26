from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 500;
kappa = 0.8;
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
const = zeros((N));
vel = zeros((N));
water = zeros((N));
for i in range(N):

    vaporData = np.load(initPath+vaporPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();
    print 'reading',vorticityPath+str(i)+'.npy';

    psi = InvPotentialVorticity(vaporData);
    v = partialX(psi);

    vaporData = (-vaporData) + 0.01;

    Q = kappa*mean(abs(v)/vaporData);
    water[i] = Q;

#    print const[i];
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    #forcing[i] =  len(flatnonzero(pvData));
    #savefig(initPath+'/figures/Forcing'+str(kappa)+'Spectrum'+str(i)+".png");
    #close();
    #figure();
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

plot(diff(water), label='kV/Q');
title('Value of $\int \kappa V/Q$  over time for k = '+str(kappa));
show();
