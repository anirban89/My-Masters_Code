from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 512;
kappa = 0.7;
initPath = '/Data/';
vorticityPath = 'EqVaportestHiRes'+str(kappa)+'/eta_data_'
vorticityPath2 = 'Vapor1.8/delta_data_'
vorticityPath3 = 'Vapor0.2/delta_data_'
vaporPath = 'Vapor2/vapor_data_'
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
forcing = zeros((512))
forcing2 = zeros((512))
forcing3 = zeros((512))

for i in range(512):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();
    pvData2 = np.load(initPath+vorticityPath2+str(i)+'.npy').transpose();
    pvData3 = np.load(initPath+vorticityPath3+str(i)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(i)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
    forcing[i] =  len(flatnonzero(pvData));
    forcing2[i] =  len(flatnonzero(pvData2));
    forcing3[i] =  len(flatnonzero(pvData3));
    time.sleep(0.05);
    #savefig(initPath+'/figures/Forcing'+str(kappa)+'Spectrum'+str(i)+".png");
    #close();
    #figure();
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

plot(forcing);
hold(True);
plot(forcing2+2);
plot(forcing3+4);
show();

#savefig(initPath+'/figures/TimeEvol'+str(kappa)+'Forcing'+str(i)+".png");
