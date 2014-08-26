from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 128;
coeff = 0.6;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
#initPath = '/home/joy/SlowCondensation/';
initPath = '/Data/';
vorticityPath = 'EqVaportestJet'+str(coeff)+'/eta_data_'
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
hold(True);
for i in arange(24,step=4):
    av = zeros((2,N));
    for j in range(5):
#pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
            pvData = np.load(initPath+vorticityPath+str(100*i+10*j)+'.npy').transpose();
            print shape(pvData);
            print 'reading',vorticityPath+str(100*i+10*j)+'.npy';
            psi = InvPotentialVorticity(pvData);
            u = -partialY(psi);
            v = partialX(psi);

            avU= average(u,0);
            av[0,:] += avU;


    print shape(av);
    print shape(x);
    fig = plot(av[0,:]/100.0,x, label='U'+str(i));
legend(loc='lower left');
figure();
for i in arange(27,step=4):
    av = zeros((2,N));
    for j in range(5):
#pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
            pvData = np.load(initPath+vorticityPath+str(100*i+10*j)+'.npy').transpose();
            print shape(pvData);
            print 'reading',vorticityPath+str(100*i+10*j)+'.npy';
            psi = InvPotentialVorticity(pvData);
            v = partialX(psi);

            avV= average(v,1);
            av[1,:] += avV;


    print shape(av);
    print shape(x);
    fig = plot(av[1,:]/100.0,x, label='V'+str(i));



#       time.sleep(0.05);
#    savefig(initPath+'EqVapor'+str(coeff)+'uVel'+str(coeff)+'-'+str(i)+".png");
#out = raw_input();
#out = int(out);

#scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
#imsave(figurePath+str(i)+".png",pvData);

legend(loc='lower left');
show();
