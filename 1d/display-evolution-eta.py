from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
from time import *;
import sys;
from scipy.stats import linregress;
import time;

N = 128;

if(len(sys.argv) < 2):
    exit;

coeff = float(sys.argv[1]);
print coeff;

if(len(sys.argv) == 3):
    time1 = int(sys.argv[2]);
else:
    time1 = 10000;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'EqVaportest512Linear'+str(coeff)+'/eta_data_'
initPath = '/home/joy/';

ion();

handle, = plot(zeros((N)));
ylim(-1,1);
draw();
hold(True);
out = 0;

offset = 00;

maxVal = zeros((time1-offset));


for i in range(offset,time1):

    j = i;

    pvData = np.load(initPath+vorticityPath+str(j)+'.npy').transpose();
    print "In time1 step: ",str(i);
#    vData = np.load(initPath+vaporPath+str(j)+'.npy').transpose();
#    deltaData = np.load(initPath+deltaPath+str(i)+'.npy').transpose();
    psi = InvPotentialVorticity(pvData);
#    rhData = (vaporField)/saturatedField;
#    vaporSeries[i] = sum(sum(vData*vData));


#    psi = InvPotentialVorticity(pvData);

#    u = partialX(psi)

    psi = psi-average(psi);
    psi = average(psi,1);
#    psi = psi/amax(abs(psi));

    #handle.set_ydata(psi);
    plot(psi);
    maxVal[j] = amax(psi);
    draw();
#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

hold(True);
figure();
slope = linregress(arange(0,len(maxVal)), log(maxVal));
plot(log(maxVal));
print slope;
title('growth rate: '+str(slope[0]));

show();


