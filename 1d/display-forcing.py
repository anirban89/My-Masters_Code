from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 512;
kappa = 0.1;
initPath = '/home/joy/EqVaportestHiRes';
vorticityPath = str(kappa)+'/delta_data_'
vaporPath = 'Vapor2/vapor_data_'
figurePath = 'Vapor2/figure_'

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
forcing = zeros((512))
#handle = imshow(zeros((N,N)),cmap=my_cmap);
#colorbar();
#draw();

time = 400;
condensing = zeros(time);
for i in range(time):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(i)+'.npy';

    condensing[i] = len(flatnonzero(pvData));
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();
#    out = fftshift(log(abs(fft2(pvData))));
#    handle.set_array(out.transpose());
#    handle.autoscale();
#    title('Spectrum with kappa='+str(kappa));
#    draw();
    #forcing[i] =  len(flatnonzero(pvData));
    #savefig(initPath+'/figures/Forcing'+str(kappa)+'Spectrum'+str(i)+".png");
    #close();
    #figure();
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

plot(condensing);
show();
