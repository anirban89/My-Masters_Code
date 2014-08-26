from coupledDerivative import *;
from boundary import *;
from finiteDifference import *;

from pylab import *;

import time;


N = 256;
i = linspace(0,N,N+2);
i = i[1:len(i)-1];

bc = boundary(yBcType='Dirichlet', bcs=( zeros((N,10)), zeros((N,10)) ));
bc1 = boundary(yBcType='Dirichlet', bcs=( zeros((N,4)), zeros((N,4)) ));

inp = exp(-pow(i-125,2)/300);
inpPrime = inp*(-2*(i-125)/300);

#inpDoublePrime
data = zeros((N,N));

for i in arange(N):
    data[i,:] = inp[:];

deriv = coupledDerivative(N,N,512.0,512.0);

for i in arange(2):
    print 'Run ', str(i);
    start = time.time();
    (outx,outx1) = deriv.partialX(data.transpose(), returnDoublePrime=True);
    boydOut = boydPartialY(data, bc1, 512.0);
    (out,out1) = deriv.partialY(data, bc, returnDoublePrime=True);
    print 'Time :', time.time()-start;

hold(True);

print shape(out);
print shape(inp);
plot(inp);
figure();
plot(inpPrime - out[1,:]);
#plot(out1[1,:]);

figure();
#plot(inpPrime);
plot(inpPrime - outx[:,1]);
figure();
plot(inpPrime - boydOut[:,1]);
#plot(outx1[:,1]);
show();
