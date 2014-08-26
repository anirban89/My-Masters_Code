from pylab import *
from spectral import *

N = 101;

x = linspace(0,2*pi, N+1);
x = x[0:-1];

print len(x);

y = linspace(0,pi,N+1);
y = y[0:-1];

print len(y)

test = zeros((N,N));

for i in arange(N):
    test[i,:] = cos(x[i])*cos(y)


out = InvPotentialVorticityCosine(test)

inverse = partialX(out,2) + partialYCosine(out,2) - out;

imshow(test.transpose())
colorbar();
figure();
imshow(inverse.transpose());
colorbar();
print norm((test-inverse),inf)
show();
