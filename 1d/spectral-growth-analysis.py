from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;

N = 391;
coeff = 0.1;
vorticityPath = '/media/FreeAgent Drive/Vapor'+str(coeff)+'/eta_data_'
vorticityPath = '/home/joy/Vapor'+str(coeff)+'/eta_data_'
#vorticityPath = '/home/joy/SlowCondensation/Vapor'+str(coeff)+'/eta_data_'
vaporPath = '/home/joymm/Vapor/vapor_data_'

k = zeros((10,N));
spectrum = zeros((N,64));

def smooth(array):
    output = zeros(len(array)-4);
    for i in arange(4,N):
        output[i-4] = sum(array[i-4:i]);

    return output;

for i in range(1,N+1):

    pvData = np.load(vorticityPath+str(i)+'.npy');
    print 'reading',vorticityPath+str(i)+'.npy';
    psi = InvPotentialVorticity(pvData);
    u = -partialY(psi);
    v = partialX(psi);
    omega = partialX(v) - partialY(u);
    omegaHat = fft(omega.transpose());
    av = average(omegaHat,axis=0);
    spectrum[i-1,:] = abs(av[0:64]);
    k[0][i-1] = abs(av[1]);
    k[1][i-1] = abs(av[2]);
    k[2][i-1] = abs(av[3]);
    k[3][i-1] = abs(av[4]);
    k[4][i-1] = abs(av[5]);
    k[5][i-1] = abs(av[20]);
    k[6][i-1] = abs(av[30]);
    k[7][i-1] = abs(av[50]);
    k[8][i-1] = abs(av[40]);
    k[9][i-1] = abs(av[12]);

index = array((1,2,3,4,5,20,30,50,40,12));
wavenumber = linspace(1,64,64);
k[k==0] = 1;
hold(True);


for i in arange(7):
#    out = smooth(k[i,:]);
    #plot(log(abs(out[i,1:N]-out[i,0:N-1])),label='wavenumber: ' + str(index[i]));
    #plot((abs(diff(k[i,:]))),label='wavenumber: ' + str(index[i]));
    plot(log((k[i,:])),label='wavenumber: ' + str(index[i]));
#legend(loc='upper right');
legend(loc='lower left');
figure();
for i in arange(4):
    t = 80*i;
    plot(log(wavenumber),log(spectrum[t,:]),label='time: '+str(t));

legend(loc='lower left');

show();

