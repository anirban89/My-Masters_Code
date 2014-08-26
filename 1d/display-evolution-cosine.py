from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 256;
initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'VaporS0.2/eta_data_'
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

handle = imshow(zeros((N/2+1,N)),cmap=my_cmap);
colorbar();
draw();
out = 1;
for i in range(1000):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(out)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(i)+'.npy';
    psi = InvPotentialVorticityCosine(pvData);
    u = -partialYCosine(psi);
    v = partialX(psi);
    omega = partialX(v) - partialY(u);
    print sqrt(sum(omega*omega));
    handle.set_array(pvData.transpose());
    handle.autoscale();
    draw();
    time.sleep(0.05);
    savefig(initPath+'CosineBeta0.2PV'+str(i)+".png");
    out = raw_input();
    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


