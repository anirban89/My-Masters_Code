from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 256;
initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'Vapor1.8/eta_data_'
vaporPath = 'Vapor2/vapor_data_'
figurePath = 'Vapor2/figure_'

ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

handle = imshow(zeros((N,N)));
colorbar();
draw();
out = 0;
for i in range(512):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();

    print shape(pvData);
    print 'reading',vorticityPath+str(i)+'.npy';
    spec = log(abs(fft2(pvData))+1);
    spec = fftshift(spec);
#    spec = fftshift(spec.transpose());
    handle.set_array(spec);
    handle.autoscale();
    draw();
    time.sleep(0.05);
    savefig(initPath+'Beta1.8PV'+str(i)+".png");

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


