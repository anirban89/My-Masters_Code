from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;

def deriveVelocities(omega1, omega2):
    '''
    omega1 = laplacian(psi1) + (psi2-psi1)
    omega2 = laplacian(psi2) + (psi1-psi2)
    '''
    sumPsi = InvLaplacian(omega1+omega2)/2.0;
    diffPsi = (omega2 - omega1)/2.0;

    psi1 = (sumPsi + diffPsi)/2.0;
    psi2 = (sumPsi - diffPsi)/2.0;

    return(psi1,psi2);


N = 512;
coeff = 0.4;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
#initPath = '/home/joy/SlowCondensation/';
initPath = '/home/joy/';
vorticityPath1 = 'TwoLayer'+str(coeff)+'/eta1_data_'
vorticityPath2 = 'TwoLayer'+str(coeff)+'/eta2_data_'

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
av = zeros((2,N));
for j in arange(0,60,5):
            pvData1= np.load(initPath+vorticityPath1+str(10*j)+'.npy').transpose();
            pvData2= np.load(initPath+vorticityPath2+str(10*j)+'.npy').transpose();
            print 'reading',vorticityPath1+str(10*j)+'.npy';
            (psi1,psi2) = deriveVelocities(pvData1,pvData2);

            u1 = -partialY(psi1);

            u2 = -partialY(psi2);

            avU= average(u1,0);
            av[0,:] += avU;

            avU= average(u2,0);
            av[1,:] += avU;



            fig = plot(av[0,:]/30.0,x, label='time: '+str(j));
            legend(loc='lower left');
show();
