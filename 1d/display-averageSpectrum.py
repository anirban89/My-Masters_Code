from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 512;
kappa = 0.9;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/Remote/';
vorticityPath = 'EqVaporHiRes'+str(kappa)+'/eta_data_'
vaporPath = 'EqVaporHiRes'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        power = power + abs(fft(field[:,i]));

    power = power/y;
    return log(power[1:N/4]);

def calcVaporPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        power = power + abs(fft(field[:,i]));

    power = power/y;
    return log(power[1:N/2]);


# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
forcing = zeros((512));
wavenumber = arange(1,N/4);
wavenumberVapor = arange(1,N/4);
hold(True);
M = 30;
out = zeros(N/4-1);
for i in range(2*M/3,M):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    value = 10*i+1;
    pvData = np.load(initPath+vorticityPath+str(value)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    out = out + calcPower(pvData);

    #draw();
    forcing[i] =  len(flatnonzero(pvData));
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

plot(log(wavenumber),out/M,label='PV');

xlabel('Wavenumber');
ylabel('Absolute Value');
title('Spectrum with kappa='+str(kappa));
#legend();
savefig(initPath+'/PV'+str(kappa)+'Spectrum.png');

out = zeros(N/4-1);
for i in range(2*M/3,M):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    value = 10*i+1;
    pvData = np.load(initPath+vaporPath+str(value)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vaporPath+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    out = out + calcPower(pvData);
    #draw();
    forcing[i] =  len(flatnonzero(pvData));
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

plot(log(wavenumberVapor),out/M,label='Vapor');
#legend();
savefig(initPath+'/Vapor'+str(kappa)+'Spectrum.png');
show();
