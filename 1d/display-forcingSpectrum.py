from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 1024;
kappa = 0.7;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/';
#initPath = '/Data/Remote/';
vorticityPath = 'EqVaporHiRes'+str(kappa)+'/eta_data_'
vaporPath = 'EqVaporHiRes'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        power = power + abs(fft(field[:,i]));

    power = power/y;
    return log(power[1:N/3]);

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
forcing = zeros((N));
wavenumber = arange(1,N/3);
wavenumberVapor = arange(1,N/2);
hold(True);
M = 5;
for i in range(M):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    value = 100*i+1;
    pvData = np.load(initPath+vorticityPath+str(value)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    out = calcPower(pvData);
    plot(log(wavenumber),out,label='time='+str(value));
    #draw();
    forcing[i] =  len(flatnonzero(pvData));
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

xlabel('Wavenumber');
ylabel('Absolute Value');
title('Spectrum with kappa='+str(kappa));
legend();
savefig(initPath+'/PVVVlarge'+str(kappa)+'Spectrum.png');
show();

for i in range(M):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    value = 100*i+1;
    pvData = np.load(initPath+vaporPath+str(value)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vaporPath+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    out = calcVaporPower(pvData);
    plot(log(wavenumberVapor),out,label='time='+str(value));
    #draw();
    forcing[i] =  len(flatnonzero(pvData));
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

xlabel('Wavenumber');
ylabel('Absolute Value');
title('Vapor Spectrum with kappa='+str(kappa));
#legend();
savefig(initPath+'/VaporVVlarge'+str(kappa)+'Spectrum.png');
show();

