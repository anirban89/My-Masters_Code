from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

N = 128;
coeff = 0.9;
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
#vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
vaporPath = 'EqVaportestHiRes'+str(coeff)+'/vapor_data_'
#initPath = '/media/FreeAgent Drive/Vapor'+str(coeff);
initPath = '/home/joy/';
#initPath = '/Data/';
figurePath = 'Figure'+str(coeff)+'/figure_'

#draw();
out = 0;

q0 = 7;
dy = 2*pi/N;
minQ = q0 - coeff*dy*N/2;

time = 3000;
U = zeros(time);
V = zeros(time);
Q = zeros(time);
P = zeros(time);

#ratioV = zeros(time);
ion()
hold(True);
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

    P[i] = sum(sum(psi*psi));
    V[i] = sum(sum(v*v));
    Q[i] = sum(sum(vaporData*vaporData));
    U[i] = sum(sum(u*u));



KE = (U+V)/amax(U+V);
V = V/amax(V);
U = U/amax(U);
Q = Q/amax(Q);
P = P/amax(P);
scatter(V,Q);
title('Scatter of V vs Q for kappa = '+str(coeff));
figure();
scatter(U,Q);
title('Scatter of U vs Q for kappa = '+str(coeff));
figure();
scatter(U,V);
title('Scatter of U vs V for kappa = '+str(coeff));
figure();
scatter(KE,Q);
title('Scatter of KE vs Q for kappa = '+str(coeff));
figure();
scatter(log(P),log(Q));
title('Scatter of log(PV) vs log(Q) for kappa = '+str(coeff));
#plot(ratioV, label = 'v^2/KE');
show();
