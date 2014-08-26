from coupledDerivative import *;
from boundary import *;
from finiteDifference import *;

from pylab import *;

import time;


N = 512;
ind = linspace(0,N,N+2);
ind = ind[1:len(ind)-1];

bc = boundary(yBcType='Dirichlet', bcs=( zeros((N,10)), zeros((N,10)) ));
bc1 = boundary(yBcType='Dirichlet', bcs=( zeros((N,4)), zeros((N,4)) ));

inpPrime = zeros((N,N));

#inpDoublePrime
data = zeros((N,N));

for i in arange(N):
        data[:,i] = exp(-pow(i-256,2)/30)*exp(-pow(ind-256,2)/30);#*exp(-pow(j-256,2)/30);
        inpPrime[:,i] = data[:,i]*(-2*(ind-256)/30);

deriv = coupledDerivative(N,N,512.0,512.0);

for i in arange(2):
    print 'Run ', str(i);
    start = time.time();
    (outx,outx1) = deriv.partialX(data, returnDoublePrime=True);
    boydOut = boydPartialX(data, 512.0);
    boydOutY = boydPartialY(data, bc1, 512.0);
    (out,out1) = deriv.partialY(data, bc, returnDoublePrime=True);
    print 'Time :', time.time()-start;


imshow(data.transpose());
colorbar();
figure();
imshow((inpPrime-outx).transpose());
colorbar();
figure();
#plot(out1[1,:]);

#plot(inpPrime);
imshow((inpPrime-boydOut).transpose());
colorbar();
'''
figure();
hold(True);
plot(outx);
figure();
plot(boydOut);
figure();
plot(inpPrime);
'''

#plot(outx1[:,1]);
show();
