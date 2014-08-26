from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 256;
coeff = 0.6;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
#initPath = '/home/joy/SlowCondensation/';
initPath = '/home/joy/';
vorticityPath = 'EqVapor'+str(coeff)+'/eta_data_'
vaporPath = 'Vapor1/vapor_data_'
figurePath = 'Vapor1/figure_'

#ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

out = 0;
x = arange(N);
for i in range(15):
    av = zeros(N);
    for j in range(10):
#pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
            pvData = np.load(initPath+vorticityPath+str(100*i+j)+'.npy').transpose();
            print shape(pvData);
            print 'reading',vorticityPath+str(100*i+j)+'.npy';
            psi = InvPotentialVorticity(pvData);
            u = -partialY(psi);
            v = partialX(psi);
            omega = partialX(v) - partialY(u);

            avPsi = average(pvData,0);
            av += avPsi;

    figure();

    print shape(av);
    print shape(x);
    fig = plot(av/100.0,x);


#       time.sleep(0.05);
    savefig(initPath+'EqVapor'+str(coeff)+'PV'+str(coeff)+'-'+str(i)+".png");
#out = raw_input();
#out = int(out);

#scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
#imsave(figurePath+str(i)+".png",pvData);


