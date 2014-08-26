from pylab import *
from spectral import *

N = 100;

x = linspace(0,2*pi, N+1);
x = x[0:-1];

print len(x);

y = linspace(0,pi,N+2);
y = y[1:-1];

print len(y)

test = zeros((N,N));

for i in arange(N):
    test[i,:] = cos(x[i])*sin(y)


out = InvPotentialVorticitySine(test)

inverse = partialX(out,2) + partialYSine(out,2) - out;

imshow(test.transpose())
colorbar();
figure();
imshow(inverse.transpose());
colorbar();
print norm((test-inverse),inf)
show();
