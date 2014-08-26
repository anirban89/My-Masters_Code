from spectral import *;
from pylab import *;

N = 1280;
t = linspace(0,8*pi,N+1);
t = t[0:len(t)-1];

test = zeros((N,N));
print 'sin', shape(sin(t)), 'test', shape(test);

for i in arange(N):
    test[i,:] = 0.2;



omega = partialY(test,order=2,length=2*pi) + \
    partialX(test,order=2, length=2*pi) - test;

new = InvPotentialVorticity(omega)

plot(test);
figure();
plot(new);
figure();
plot(new-test);
print test
print new
show();


