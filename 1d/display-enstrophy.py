from pylab import *
from boundary import *
from finiteDifference import *

N = 15;
dim = 350;
north = zeros((dim,4));
south = zeros((dim,4));
bc = boundary(yBcType='Dirichlet',bcs=(north, south));
hold(True);
ion();
y = zeros(dim);
for i in arange(N):
    y[i] = (i-N/2+1)*2*pi/dim;

epsInv = 10;
numerator = zeros((dim,dim));
denominator = zeros((dim,dim));

for i in arange(1,N):
        u = \
        np.load('/home/joy/EqBeta/uVel_data_'+str(i)+'.npy').transpose();
        v = \
        np.load('/home/joy/EqBeta/vVel_data_'+str(i)+'.npy').transpose();
        n = \
        np.load('/home/joy/EqBeta/eta_data_'+str(i)+'.npy').transpose();

        plot(average(n[170:180,:], axis=1));
        omega = boydPartialX(v) - boydPartialY(u,bc);
        for i in arange(dim):
                numerator[:,i] = epsInv*y[i] + omega[:,i];
        denominator = epsInv + n;


        print str(i)+': enstrophy : ', sum(omega*omega), \
        'energy : ', sum(u*u + v*v + 10*n*n), \
        'PV:', sum(numerator/denominator);

show();
