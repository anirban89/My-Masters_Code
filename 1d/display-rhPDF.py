from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
from time import *;

N = 128;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
coeff = 0.9;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/';
#initPath = '/Data/';
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
vaporPath = 'EqVaportestHiRes'+str(coeff)+'/vapor_data_'

ion()
# making my own colormap

x = array([-.5,-.4,-.3,-.2,-.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]);

handle, = plot(x,2*x,'o');
time = 3093;
for i in range(time):
    out = 10*i;

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    vaporData = np.load(initPath+vaporPath+str(out)+'.npy');
    print 'reading',vaporPath+str(out)+'.npy';
    
#    (n,bins) = histogram(vaporData,bins=[-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05],normed=True);
    (n,bins) = \
    histogram(vaporData,bins=[-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,1.1],normed=True);
    print shape(n), shape(bins);
    handle.set_ydata(n);
    draw();


#    savefig(initPath+'vaporCosinePDF'+str(i)+".png");

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);


