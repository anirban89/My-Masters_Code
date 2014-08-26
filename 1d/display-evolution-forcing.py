from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 512;
coeff = 0.7;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'EqVaporHiResNoForcing'+str(coeff)+'/eta_data_'
#initPath = '/home/joy/';
initPath = '/Data/';
#vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
vaporPath = 'EqVaportestHiRes'+str(coeff)+'/vapor_data_'
deltaPath = 'EqVaportestHiRes'+str(coeff)+'/delta_data_'
figurePath = 'Figure'+str(coeff)+'/figure_'

q0 = 7;
dy = 2*pi/N;

saturatedField = zeros((N,N));

for i in arange(N/2):
    saturatedField[:,i] = q0 - coeff*dy*(N/2-i);
    saturatedField[:,N-i-1] = q0 - coeff*dy*(N/2-i);

ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

handle1 = imshow(zeros((N,N)));
colorbar();
draw();
out = 0;
time = 132;

vaporSeries = zeros(time);

for i in range(time):

    j = 10*i;

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    vData = np.load(initPath+vaporPath+str(j)+'.npy').transpose();
    deltaData = np.load(initPath+deltaPath+str(j)+'.npy').transpose();
    print 'reading',deltaPath+str(j)+'.npy';
    vData = vData*deltaData;
    handle1.set_array(vData.transpose());
#    vaporSeries[i] = len(flatnonzero(vData));
    vaporSeries[i] = sum(sum(vData));
    handle1.autoscale();
    draw();
   # time.sleep(0.05);
#    savefig(initPath+figurePath+'Forcing'+str(coeff)+'-'+str(i)+".png");
#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

figure();
plot(vaporSeries);
show();
