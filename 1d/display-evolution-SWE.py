from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 256;
coeff = 0.6;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'SWE'+str(coeff)+'/eta_data_'
#initPath = '/media/FreeAgent Drive/Vapor'+str(coeff);
initPath = '/home/joy/';
vaporPath = 'Vapor2/vapor_data_'
figurePath = 'VaporFigure'+str(coeff)+'/figure_'

ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

handle = imshow(zeros((N,N)),cmap=my_cmap);
colorbar();
#draw();
out = 1;
for i in range(600):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(out)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(out)+'.npy';
    psi = partialX(pvData, 2) + partialY(pvData, 2) - pvData;
    handle.set_array(psi.transpose());
    handle.autoscale();
    draw();
   # time.sleep(0.05);
   # savefig(initPath+figurePath+'PV'+str(coeff)+'-'+str(i)+".eps");
    out = raw_input();
    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


