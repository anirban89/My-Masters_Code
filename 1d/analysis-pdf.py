from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;

vorticityPath = '/home/joymm/Vapor/eta_data_'
vaporPath = '/home/joymm/Vapor/vapor_data_'
N = 290

pvData = np.load(vorticityPath+str(N)+'.npy');
print 'reading',vorticityPath+str(N)+'.npy';
psi = InvPotentialVorticity(pvData);
u = -partialY(psi);
v = partialX(psi);
omega = partialX(v) - partialY(u);
variance = var(omega);
(n,bins) = histogram(omega,bins=50,normed=True);
pdf = zeros((N-1,len(n)));
bins = zeros((N-1,len(bins)));
print bins

for i in range(1,N):

    pvData = np.load(vorticityPath+str(i)+'.npy');
    print 'reading',vorticityPath+str(i)+'.npy';
    psi = InvPotentialVorticity(pvData);
    u = -partialY(psi);
    v = partialX(psi);
    omega = partialX(v) - partialY(u);
    variance = var(omega);
    (pdf[i-1,:],bins[i-1,:]) = histogram(omega,bins=50,normed=True);



hold(True);
for i in arange(N):
    plot(.5*(bins[i,1:]+bins[i,:-1]),pdf[i,:],label='time: ' + str(i));
    legend(loc='lower left');
    savefig('/home/joymm/Vapor/figure'+str(i+1)+'.png');
    print 'savefig(/home/joymm/Vapor/figure'+str(i)+'.png';
    figure();

#show();

