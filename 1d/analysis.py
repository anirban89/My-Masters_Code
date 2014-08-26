from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;

vorticityPath = '/home/joymm/Vapor/vorticity_data_'
vaporPath = '/home/joymm/Vapor/vapor_data_'


energy = zeros(121);

for i in range(1,121):

    pvData = np.load(vorticityPath+str(i)+'.npy');
    print 'reading',vorticityPath+str(i)+'.npy';
    psi = InvPotentialVorticity(pvData);
    u = -partialY(psi);
    v = partialX(psi);
    kinetic = psi*(u*u + v*v);
    potential = 9.8*psi*psi;
    energy[i] = sum(kinetic + potential);


plot(log(energy));
show();

