from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 512;
coeff = 0.7;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
#vorticityPath = 'EqVaporSuperHiRes'+str(coeff)+'/eta_data_'
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
#initPath = '/media/FreeAgent Drive/Vapor'+str(coeff);
#initPath = '/home/joy/Remote/';
initPath = '/home/joy/';
#initPath = '/Data/';
figurePath = 'Figure'+str(coeff)+'/figure_'

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
#draw();
out = 0;
time = 196;
for i in range(time):

    j = 5*i;
    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(j)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(j)+'.npy';
    psi = InvPotentialVorticity(pvData);
    u = -partialY(psi)
    handle.set_array(u.transpose());
    handle.autoscale();
    draw();
#    time.sleep(0.05);
#    savefig(initPath+figurePath+'PV'+str(coeff)+'-'+str(i)+".png");
#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


