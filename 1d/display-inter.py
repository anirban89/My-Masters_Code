from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 256;
vorticityPath = '/home/joymm/VaporS4/inter_data_'
vaporPath = '/home/joymm/VaporS1/vapor_data_'
figurePath = '/home/joymm/VaporS1/figure_'

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
out = 2;
for i in range(1000):

    #pvData = np.load(vorticityPath+str(out)+'.npy').transpose();
    pvData = np.load(vorticityPath+str(out)+'.npy');
    print shape(pvData);
    print 'reading',vorticityPath+str(out)+'.npy';
    #psi = InvPotentialVorticityCosine(pvData);
    #u = -partialYCosine(psi);
    #v = partialX(psi);
    #omega = partialX(v) - partialY(u);
    #print sqrt(sum(pvData*pvData));
    handle.set_array(pvData.transpose());
    handle.autoscale();
    draw();
    time.sleep(0.1);
    out = raw_input();
    out = int(out);
#    savefig(figurePath+str(i)+".png");

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


