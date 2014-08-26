from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;
import sys;


if (len(sys.argv) < 2):
        print 'please give value of kappa as argument';
        exit();

kappa = float(sys.argv[1]);
print kappa;

if (len(sys.argv) == 3):
        M = int(sys.argv[2]);
else:
        M = 10000;

N = 128;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/';
#initPath = '/Data/';
vorticityPath = 'EqVaportestHiRes'+str(kappa)+'/eta_data_'
vaporPath = 'EqVaporHiRes'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

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
hold(True);
forcing = zeros((M,N/3));
for i in range(M):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    value = 10*i;
    pvData = np.load(initPath+vorticityPath+str(value)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    out = calcPower(pvData);
    #draw();
    forcing[i,:] =  out;
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

figure(figsize=(10,8));
plot(forcing[:,0], label='WN=0');
plot(forcing[:,1], label='WN=1');
plot(forcing[:,2], label='WN=2');
plot(forcing[:,3], label='WN=3');
plot(forcing[:,4], label='WN=4');
title('Spectrum with kappa='+str(kappa));
legend(loc='lower right');

#legend();
savefig(initPath+'/ForcingSpectrumZonal'+str(kappa)+'.png');
#show();

