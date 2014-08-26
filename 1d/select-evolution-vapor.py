from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import sys;

N = 128;

if(len(sys.argv) < 2):
    exit;

coeff = float(sys.argv[1]);
print coeff;

#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
initPath = '/home/joy/';
#initPath = '/Data/';
#vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
#vaporPath = 'EqVaportestHiResBoundaryForcing'+str(coeff)+'/vapor_data_'
#vaporPath = 'EqVaportestLoRes'+str(coeff)+'/vapor_data_'
vaporPath = 'EqVaportestHiRes'+str(coeff)+'/vapor_data_'
deltaPath = 'Vapor'+str(coeff)+'/delta_data_'
figurePath = 'Figure'+str(coeff)+'/figure_'

q0 = 7;
dy = 2*pi/N;

saturatedFieldG = zeros((N,N));

for i in arange(N/2):
    saturatedFieldG[:,i] = q0 - coeff*dy*(N/2-i);
    saturatedFieldG[:,N-i-1] = q0 - coeff*dy*(N/2-i);

# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

out = 0;

offset = 000;

for i in arange (3650,3850):

    print 'Enter value of time ';
#    i = int(input());
    j = 10*i;

    pvData = np.load(initPath+vorticityPath+str(j)+'.npy').transpose();
#    deltaData = np.load(initPath+deltaPath+str(i)+'.npy').transpose();
    print 'reading',vaporPath+str(j)+'.npy';
#    rhData = (vaporField)/saturatedField;
#    vaporSeries[i] = sum(sum(vData*vData));


    psi = InvPotentialVorticity(pvData);

    u = -partialY(psi);
    v = partialX(psi);
    
    figure();
    imshow(u.transpose());
    colorbar();
    savefig(initPath+'/figs/'+'uVel'+str(j)+'.png');
    figure();
    imshow(v.transpose());
    colorbar();
    savefig(initPath+'/figs/'+'vVel'+str(j)+'.png');

#    show();

    print 'Save? '
 #   answer = raw_input();

#    if(answer == 'y'):

    print 'Saving figure';
#    imsave(initPath+'/figs/'+str(coeff)+'uVel'+str(j)+'.png',u.transpose());
#    imsave(initPath+'/figs/'+str(coeff)+'vVel'+str(j)+'.png',v.transpose());


   # time.sleep(0.05);
#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

exit();

