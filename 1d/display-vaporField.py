from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 256;
#initPath = '/media/FreeAgent Drive/Vapor-highres/';
coeff = 0.6;
initPath = '/home/joy/';
vorticityPath = 'Vapor2/eta_data_'
vaporPath = 'EqVapor'+str(coeff)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

ion()
# making my own colormap
cdict_vapor = {
    'red'  :  ((0., 0, 0), (0.5, 0.5, 0.5), (1., 0.5, 1)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 0.5, 0), (0.5, 0.5, 0.5), (1, 0.5, 1))
}
my_cmap_vapor = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict_vapor, 1024)

handle = imshow(zeros((N,N)));
colorbar();
draw();


out = 0;
x = arange(512);
for i in range(1000):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    vaporData = np.load(initPath+vaporPath+str(out)+'.npy');
    print 'reading',vaporPath+str(out)+'.npy';
    
    handle.set_array(vaporData);
    handle.autoscale();
    draw();

#    time.sleep(0.05);
#    savefig(initPath+'vaporLoBeta'+str(i)+".png");
    out = raw_input();
    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


