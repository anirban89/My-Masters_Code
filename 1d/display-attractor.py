from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import axes3d, Axes3D;

N = 128;
coeff = 0.9;
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
#vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
vaporPath = 'EqVaportestHiRes'+str(coeff)+'/vapor_data_'
#initPath = '/media/FreeAgent Drive/Vapor'+str(coeff);
initPath = '/home/joy/';
#initPath = '/Data/';
figurePath = 'Figure'+str(coeff)+'/figure_'

fig = plt.figure();
ax = Axes3D(fig);
#draw();
out = 0;

q0 = 7;
dy = 2*pi/N;
minQ = q0 - coeff*dy*N/2;

time = 1413;
U = zeros(time);
V = zeros(time);
Q = zeros(time);
#ratioV = zeros(time);
for i in range(time):

    j = 10*i;

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(j)+'.npy').transpose();
    vaporData = np.load(initPath+vaporPath+str(j)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(j)+'.npy';
    psi = InvPotentialVorticity(pvData);
    v = partialX(psi);
    u = -partialY(psi);

    Q[i] = sum(sum(vaporData*vaporData));

    U[i] = sum(sum(real(u*u)));
    V[i] = sum(sum(real(v*v)));
#    ratioV[i] = sum(sum(real(v*v)))/sum(sum(real(u*u + v*v)));
    #ratio[i] = sum(sum(real(10*psi*psi)))/totalPotentialEnergy;

#    ratio[i] = sum(sum(real(u*u + v*v)));

#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

ax.plot(U,Q, 'o');
title('Attractor' );
#plot(ratioV, label = 'v^2/KE');
show();
